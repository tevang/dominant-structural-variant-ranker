"""fix-auto3d-integration §5: consolidated Auto3D failure reporting.

- One root-cause record per distinct failure; affected candidates carry
  short references (task 5.1).
- Terminal failures are remembered per stage and skip subsequent identical
  invocations (task 5.2).
- Per-candidate warning text is bounded (task 5.3).
"""

import csv
import json
from pathlib import Path

from rdkit import Chem

from dsvr.chemistry import tautomer_auto3d_filter as tautomer_filter
from dsvr.chemistry.tautomer_auto3d_filter import _Candidate, filter_tautomers_with_auto3d
from dsvr.config import RunConfig
from dsvr.models import ProtomerRecord, make_tautomer_id
from dsvr.reporting.auto3d_diagnostics import _BOOKS, bounded, failure_book_for
from dsvr.runners.auto3d_runner import Auto3DExecutionError

SEMLOCK_ERROR = (
    "RuntimeError: A SemLock created in a fork context is being shared "
    "with a process in a spawn context."
)


def _protomer(protomer_id: str, smiles: str = "CC(=O)C") -> ProtomerRecord:
    mol = Chem.MolFromSmiles(smiles)
    return ProtomerRecord(
        id=protomer_id,
        parent_id=f"{protomer_id}_parent",
        input_molecule_id=f"{protomer_id}_input",
        molname=protomer_id,
        canonical_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        molecular_formula="C3H6O",
        formal_charge=0,
        explicit_proton_count=6,
        source_software="test",
        protomer_index=1,
        rdkit_mol=mol,
    )


def _candidate(protomer: ProtomerRecord, index: int, smiles: str) -> _Candidate:
    mol = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    metadata = {"auto3d_tautomer_filtering": True, "tautomer_filtering_stage": "pre_stereo"}
    return _Candidate(
        index=index,
        tautomer_id=make_tautomer_id(protomer.id, index, canonical, isomeric, metadata),
        molecule=mol,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        rdkit_score=0.0,
        is_input_tautomer=index == 1,
        is_canonical_tautomer=index == 1,
    )


def _two_candidates(protomer: ProtomerRecord) -> list[_Candidate]:
    return [
        _candidate(protomer, 1, "CC(=O)C"),
        _candidate(protomer, 2, "C=C(O)C"),
    ]


def test_distinct_failure_written_once_with_short_candidate_refs(tmp_path, monkeypatch) -> None:
    """Task 5.1: a multi-protomer stage failing identically produces exactly
    one root-cause record; candidate rows reference it with short notes."""

    _BOOKS.clear()
    invocations: list[str] = []

    def fail_auto3d(input_path: Path, output_dir: Path, **kwargs):
        invocations.append(kwargs["model"])
        raise Auto3DExecutionError("mock global outage: auto3d exploded identically")

    monkeypatch.setattr(tautomer_filter, "run_auto3d", fail_auto3d)
    protomers = [_protomer("mol_a"), _protomer("mol_b"), _protomer("mol_c")]
    candidates_by_id = {p.id: _two_candidates(p) for p in protomers}
    monkeypatch.setattr(
        tautomer_filter,
        "_enumerate_candidates",
        lambda protomer, config: (candidates_by_id[protomer.id], None),
    )
    config = RunConfig(output_dir=tmp_path / "run")

    filter_tautomers_with_auto3d(protomers, config)

    jsonl = tmp_path / "run" / "auto3d_root_causes.jsonl"
    records = [
        json.loads(line)
        for line in jsonl.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    # One distinct root cause; occurrences are appended with running counts.
    distinct = {record["root_cause_id"] for record in records}
    assert len(distinct) == 1
    assert records[-1]["occurrences"] > 1

    selected_csv = tmp_path / "run" / "enumeration" / "tautomers" / "tautomers_selected.csv"
    rows = list(csv.DictReader(selected_csv.open()))
    assert rows
    refs = set()
    for row in rows:
        warnings_text = row["warnings"]
        assert len(warnings_text) <= 3 * 300, warnings_text[:400]
        assert "auto3d_failed:EXECUTION_ERROR" in warnings_text
        assert "ref " in warnings_text
        refs.add(warnings_text.split("ref ")[1].rstrip("() ")[:12])
    assert len(refs) == 1, refs


def test_terminal_failure_memory_skips_identical_invocations(tmp_path, monkeypatch) -> None:
    """Task 5.2/5.3: after one SemLock crash the stage remembers the terminal
    infra failure — remaining protomers use the fallback directly and no
    second identical invocation happens."""

    _BOOKS.clear()
    invocations: list[str] = []

    def semlock_auto3d(input_path: Path, output_dir: Path, **kwargs):
        invocations.append(kwargs["model"])
        raise Auto3DExecutionError(SEMLOCK_ERROR)

    monkeypatch.setattr(tautomer_filter, "run_auto3d", semlock_auto3d)
    protomers = [_protomer(f"mol_{index}") for index in range(4)]
    candidates_by_id = {p.id: _two_candidates(p) for p in protomers}
    monkeypatch.setattr(
        tautomer_filter,
        "_enumerate_candidates",
        lambda protomer, config: (candidates_by_id[protomer.id], None),
    )
    config = RunConfig(output_dir=tmp_path / "run")

    filter_tautomers_with_auto3d(protomers, config)

    assert invocations == ["ANI2xt"], invocations
    jsonl = tmp_path / "run" / "auto3d_root_causes.jsonl"
    records = [
        json.loads(line)
        for line in jsonl.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert len(records) == 1
    assert records[0]["error_class"] == "INFRA_MULTIPROCESSING"
    stage_note = (
        tmp_path / "run" / "enumeration" / "tautomers" / "tautomers_selected.csv"
    ).read_text(encoding="utf-8")
    assert "skipped per failure memory" not in stage_note  # only attempt rows carry it
    assert "auto3d_failed:INFRA_MULTIPROCESSING" in stage_note


def test_engine_incompatible_memory_is_engine_scoped(tmp_path) -> None:
    _BOOKS.clear()
    book = failure_book_for(tmp_path / "run")
    first = book.record_failure("stage", "ANI2xt", "Only AIMNET can handle: ['x']")
    assert first.error_class == "ENGINE_INCOMPATIBLE"
    assert book.terminal_reference("stage", "ANI2xt") is not None
    assert book.terminal_reference("stage", "AIMNET") is None


def test_bounded_note_truncates_long_text() -> None:
    assert bounded("short") == "short"
    long_text = "x" * 1000
    result = bounded(long_text, max_chars=300)
    assert len(result) == 300
    assert result.endswith("…")


def test_failures_differing_only_in_uuid_and_paths_dedup_to_one_cause(tmp_path) -> None:
    """Escalation-review regression: normalization must cover every line of
    the excerpt, not just the first."""

    _BOOKS.clear()
    book = failure_book_for(tmp_path / "run")
    base = (
        "Auto3D failed. Tried commands:\n"
        "python /tmp/{pid}_auto3d_tautomers_{uid}/auto3d_ANI2xt/_auto3d_v3_wrapper.py run "
        "--job-name final_{uid} exited 1: boom"
    )
    first = book.record_failure(
        "stage", "ANI2xt", base.format(pid="mol_a", uid="aaaa1111")
    )
    second = book.record_failure(
        "stage", "ANI2xt", base.format(pid="mol_b", uid="bbbb2222")
    )
    assert first.root_cause_id == second.root_cause_id
    assert second.count == 2
    lines = [
        json.loads(line)
        for line in (tmp_path / "run" / "auto3d_root_causes.jsonl").read_text().splitlines()
    ]
    assert lines[-1]["occurrences"] == 2
    open_ids = {line["root_cause_id"] for line in lines}
    assert len(open_ids) == 1
