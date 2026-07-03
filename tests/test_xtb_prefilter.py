from __future__ import annotations

from pathlib import Path

from dsvr.config import RunConfig
from dsvr.filtering import xtb_prefilter
from dsvr.filtering.xtb_prefilter import apply_xtb_prefilter, write_xtb_prefilter_outputs
from dsvr.models import SeedConformerRecord
from dsvr.runners.xtb_prefilter_runner import XtbPrefilterResult


def test_xtb_prefilter_can_be_disabled() -> None:
    seeds = [_seed(1, 0.0), _seed(2, 1.0)]

    selected, decisions = apply_xtb_prefilter(seeds, RunConfig())

    assert selected == seeds
    assert all(decision.selected for decision in decisions)
    assert {decision.reason for decision in decisions} == {"xtb_prefilter_disabled"}


def test_xtb_prefilter_prunes_high_energy_variants(tmp_path: Path, monkeypatch) -> None:
    seeds = [_seed(1, 0.0), _seed(2, 1.0), _seed(3, 2.0)]
    energies = {seeds[0].id: -100.0, seeds[1].id: -95.0, seeds[2].id: -80.0}

    def fake_run(seed, config):
        return XtbPrefilterResult(
            seed_id=seed.id,
            parent_stereo_id=seed.parent_id,
            input_molecule_id=seed.input_molecule_id,
            energy_kcal_mol=energies[seed.id],
            workdir=tmp_path / seed.id,
            warnings=[],
        )

    monkeypatch.setattr(xtb_prefilter, "run_xtb_prefilter", fake_run)
    config = RunConfig(
        output_dir=tmp_path,
        xtb_prefilter={
            "enabled": True,
            "keep_within_kcal_mol": 10.0,
            "keep_top_n_per_molecule": 2,
            "max_variants_per_molecule": 2,
        },
    )

    selected, decisions = apply_xtb_prefilter(seeds, config)
    write_xtb_prefilter_outputs(tmp_path / "filtering", decisions)

    assert len({seed.parent_id for seed in selected}) == 2
    assert any(not decision.selected for decision in decisions)
    assert (tmp_path / "filtering" / "xtb_prefilter_decisions.csv").exists()


def _seed(index: int, energy: float) -> SeedConformerRecord:
    return SeedConformerRecord(
        id=f"seed_{index}",
        parent_id=f"stereo_{index}",
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles="CCO",
        isomeric_smiles="CCO",
        molecular_formula="C2H6O",
        formal_charge=0,
        explicit_proton_count=6,
        source_software="test",
        source_python_function="test",
        conformer_index=1,
        energy_kcal_mol=energy,
    )
