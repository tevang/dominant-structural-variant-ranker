"""Regression tests for shared microstates across molecules in one Uni-Pka batch.

Uni-Pka 0.3.2 bug (fixed in the vendored ``containers/unipka.py``): the streaming pipeline
identified each microstate by its plain SMILES and attributed it to a single parent molecule
(``microstate_to_smi[ms] = smi`` — last writer wins). When two input molecules of one batch
shared microstates by exact SMILES collision (e.g. the conjugate pair aniline.[Cl-] /
anilinium.[Cl-], which share both the neutral form and the cation), the shared microstates
were credited twice to the last-enumerated molecule and zero times to the first. The first
molecule never reached its microstate count, was silently abandoned in ``pending``, and the
run exited 0 with fewer output lines than input lines — downstream tools then kept the
unprotonated input state without any warning.

The tests below drive the vendored ``UnipkaStream`` directly with a fake GPU predictor and a
sequential pool, so they need neither the container nor the model weights.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
# overridable only to let maintenance scripts run these tests against a modified copy
VENDORED_SCRIPT = Path(os.environ.get("DSVR_UNIPKA_SCRIPT", REPO_ROOT / "containers" / "unipka.py"))

ACID_SMI = "c1ccccc1N.[Cl-]"     # aniline hydrochloride (neutral base form)
BASE_SMI = "c1ccccc1[NH+].[Cl-]"  # anilinium chloride (cationic form)
NEUTRAL_MS = "c1ccccc1N.[Cl-]"
CATION_MS = "c1ccccc1[NH+].[Cl-]"


def _load_vendored_unipka() -> ModuleType | None:
    try:
        import torch  # noqa: F401  (imported at module level by the vendored script)
    except ImportError:
        return None
    if "vendored_unipka" in sys.modules:
        return sys.modules["vendored_unipka"]
    spec = importlib.util.spec_from_file_location("vendored_unipka", VENDORED_SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["vendored_unipka"] = module
    spec.loader.exec_module(module)
    return module


unipka = pytest.importorskip("torch", reason="vendored unipka.py imports torch")
vendored_unipka = _load_vendored_unipka()
if vendored_unipka is None:  # pragma: no cover - defensive
    pytest.skip("could not load vendored containers/unipka.py", allow_module_level=True)


class _SequentialPool:
    """Stand-in for multiprocessing.Pool that executes tasks in-process."""

    def __init__(self, *_args, **_kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        return False

    @staticmethod
    def imap_unordered(fn, iterable, chunksize=1):
        return [fn(item) for item in iterable]


class _FakePredictor:
    """Deterministic stand-in for FreeEnergyPredictor keyed by microstate SMILES."""

    def __init__(self, energies: dict[str, float]):
        self._energies = energies
        self.predicted: list[str] = []

    def predict(self, smiles_list):
        self.predicted.extend(smiles_list)
        return {smi: self._energies[smi] for smi in smiles_list if smi in self._energies}


def _hand_crafted_ensemble(smi: str, template_a2b, template_b2a):
    """Two-molecule toy batch in which both molecules share both microstates exactly."""
    if smi == ACID_SMI:
        ensemble = {0: [NEUTRAL_MS], 1: [CATION_MS]}
    elif smi == BASE_SMI:
        ensemble = {1: [CATION_MS], 0: [NEUTRAL_MS]}
    else:  # pragma: no cover - the tests only feed the pair above
        raise AssertionError(f"unexpected SMILES {smi}")
    return smi, ensemble


@pytest.fixture()
def stream_env(monkeypatch):
    monkeypatch.setattr(vendored_unipka, "Pool", _SequentialPool)
    monkeypatch.setattr(vendored_unipka, "get_ensemble", _hand_crafted_ensemble)


def _run_stream(predictor, source, gpu_trigger_microstates=10_000):
    pipeline = vendored_unipka.UnipkaStream(
        template_a2b=None,
        template_b2a=None,
        predictor=predictor,
        patterns=[],
        ncpu=1,
        pH=7.4,
        n_forms=1,
        min_occupancy=0.0,
        # default: single flush at the end, as in a small batch
        gpu_trigger_microstates=gpu_trigger_microstates,
        gpu_trigger_timeout=3600.0,
    )
    return list(pipeline.process(source))


def test_conjugate_pair_in_one_batch_keeps_both_molecules(stream_env):
    """Both members of a pair sharing all microstates must produce a MolResult."""
    # the cation must beat the pH/charge reweighting term ln10*(pH-6.5)~2.07 to be major at pH 7.4
    predictor = _FakePredictor({NEUTRAL_MS: 0.0, CATION_MS: -5.0})
    source = iter([(ACID_SMI, "M1"), (BASE_SMI, "M2")])

    results = _run_stream(predictor, source)

    by_name = {res.name: res for res in results}
    assert set(by_name) == {"M1", "M2"}, (
        "one molecule of the conjugate pair was silently dropped from the batch output"
    )
    # each molecule must own both microstates in its ensemble
    for res in results:
        microstates = {ms for members in res.ensemble_free_energy.values() for ms, _e in members}
        assert microstates == {NEUTRAL_MS, CATION_MS}
    # the cation is the most populated form at pH 7.4 for both (fake energies above)
    assert by_name["M1"].forms[0][0] == CATION_MS
    assert by_name["M2"].forms[0][0] == CATION_MS


def test_shared_microstates_predicted_once(stream_env):
    """A microstate shared by several molecules must be sent to the GPU only once."""
    predictor = _FakePredictor({NEUTRAL_MS: 0.0, CATION_MS: -2.0})
    source = iter([(ACID_SMI, "M1"), (BASE_SMI, "M2"), (ACID_SMI, "M1b")])

    results = _run_stream(predictor, source)

    assert sorted(res.name for res in results) == ["M1", "M1b", "M2"]
    assert sorted(predictor.predicted) == sorted([NEUTRAL_MS, CATION_MS])


def test_duplicate_microstate_within_one_molecule_completes(stream_env, monkeypatch):
    """A microstate listed twice within one molecule's ensemble must not block completion."""
    duplicated = {0: [NEUTRAL_MS, NEUTRAL_MS], 1: [CATION_MS]}
    monkeypatch.setattr(
        vendored_unipka,
        "get_ensemble",
        lambda smi, template_a2b, template_b2a: (smi, duplicated),
    )
    predictor = _FakePredictor({NEUTRAL_MS: 0.0, CATION_MS: -2.0})

    results = _run_stream(predictor, iter([(ACID_SMI, "M1")]))

    assert [res.name for res in results] == ["M1"]
    assert results[0].forms, "molecule with a duplicated microstate yielded no forms"


def test_late_molecule_credited_from_earlier_flush(stream_env):
    """A molecule whose shared microstates were predicted in a previous flush must still complete.

    With the GPU trigger set to fire after the first molecule, the pair is flushed separately:
    the second molecule registers only after the first (and the shared microstates' prediction)
    is already done, so it must be credited from the stored predictions rather than hang.
    """
    predictor = _FakePredictor({NEUTRAL_MS: 0.0, CATION_MS: -5.0})
    source = iter([(ACID_SMI, "M1"), (BASE_SMI, "M2")])

    results = _run_stream(predictor, source, gpu_trigger_microstates=2)

    by_name = {res.name: res for res in results}
    assert set(by_name) == {"M1", "M2"}, (
        "a molecule registered after a mid-stream flush of its shared microstates was dropped"
    )
    # each unique microstate is still predicted exactly once across both flushes
    assert sorted(predictor.predicted) == sorted([NEUTRAL_MS, CATION_MS])


def test_duplicate_smiles_redelivery_does_not_double_emit(stream_env):
    """A repeated input SMILES delivered again after completion must not re-emit its rows.

    The priority stream does not deduplicate, so the same SMILES can be delivered several
    times. With the GPU trigger firing on the first delivery, the molecule completes before
    the later deliveries arrive; each input row must still produce exactly one output row and
    no microstate may be predicted twice.
    """
    predictor = _FakePredictor({NEUTRAL_MS: 0.0, CATION_MS: -5.0})
    source = iter([(ACID_SMI, "M1"), (ACID_SMI, "M1"), (ACID_SMI, "M1")])

    results = _run_stream(predictor, source, gpu_trigger_microstates=2)

    # one output row per input row (upstream smi_to_names semantics) — not a second full batch
    assert [res.name for res in results] == ["M1", "M1", "M1"], (
        f"molecule delivered 3 times produced {len(results)} rows"
    )
    assert sorted(predictor.predicted) == sorted([NEUTRAL_MS, CATION_MS]), (
        "redelivered molecule re-predicted its microstates"
    )
