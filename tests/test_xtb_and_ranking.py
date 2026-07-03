import math
import subprocess
from pathlib import Path

import pytest

from dsvr.config import RunConfig
from dsvr.models import CrestConformerRecord, ThermoRecord
from dsvr.parsing.xtb_outputs import parse_xtb_energy, parse_xtb_thermo
from dsvr.ranking.boltzmann import R_KCAL_MOL_K
from dsvr.ranking.population import (
    CROSS_PROTOMER_WARNING,
    compute_delta_g_and_populations,
    write_ranked_outputs,
)
from dsvr.runners import xtb_runner
from dsvr.runners.xtb_runner import run_xtb_thermo
from dsvr.utils.units import hartree_to_kcal_mol


def test_parse_xtb_energy_sample_log() -> None:
    parsed = parse_xtb_energy(Path("tests/data/xtb_energy.log"))

    assert parsed.electronic_energy_hartree == -40.0
    assert parsed.electronic_energy_kcal_mol == pytest.approx(hartree_to_kcal_mol(-40.0))


def test_parse_xtb_thermo_sample_log() -> None:
    parsed = parse_xtb_thermo(Path("tests/data/xtb_thermo.log"))

    assert parsed.electronic_energy_hartree == -40.0
    assert parsed.gibbs_free_energy_hartree == -39.98
    assert parsed.gibbs_free_energy_kcal_mol == pytest.approx(hartree_to_kcal_mol(-39.98))
    assert parsed.enthalpy_hartree == -39.99
    assert parsed.entropy_cal_mol_k == 21.5


def test_boltzmann_populations_are_numerically_correct() -> None:
    config = RunConfig()
    records = [
        _thermo("a", "mol", "C2H6O", 6, 0.0),
        _thermo("b", "mol", "C2H6O", 6, 1.0),
    ]

    ranked = compute_delta_g_and_populations(records, config)

    expected_minor = math.exp(-1.0 / (R_KCAL_MOL_K * config.chemistry.temperature_kelvin))
    expected_major_pop = 1.0 / (1.0 + expected_minor)
    expected_minor_pop = expected_minor / (1.0 + expected_minor)
    assert ranked[0].relative_free_energy_kcal_mol == 0.0
    assert ranked[0].boltzmann_population == pytest.approx(expected_major_pop)
    assert ranked[1].relative_free_energy_kcal_mol == 1.0
    assert ranked[1].boltzmann_population == pytest.approx(expected_minor_pop)


def test_population_groups_same_formula_are_separate() -> None:
    config = RunConfig(thermo={"population_scope": "same_formula"})
    records = [
        _thermo("a", "mol1", "C2H6O", 6, 0.0),
        _thermo("b", "mol1", "C2H6O", 6, 1.0),
        _thermo("c", "mol2", "C2H7N", 7, 100.0),
    ]

    ranked = compute_delta_g_and_populations(records, config)

    isolated = next(record for record in ranked if record.parent_id == "c")
    assert isolated.relative_free_energy_kcal_mol == 0.0
    assert isolated.boltzmann_population == 1.0
    assert not isolated.warnings


def test_population_scope_same_charge_warns_for_mixed_formula() -> None:
    config = RunConfig(thermo={"population_scope": "same_charge"})
    records = [
        _thermo("a", "mol1", "C2H6O", 6, 0.0, charge=0),
        _thermo("b", "mol2", "C2H7N", 7, 1.0, charge=0),
    ]

    ranked = compute_delta_g_and_populations(records, config)

    assert len(ranked) == 2
    assert all(record.population_scope == "same_charge" for record in ranked)
    assert all(record.approximate_population for record in ranked)
    assert all(CROSS_PROTOMER_WARNING in record.warnings for record in ranked)
    assert all(
        record.metadata["ranking"]["mixed_formula_or_proton_count"] for record in ranked
    )


def test_population_warning_for_mixed_formula_all_approximate() -> None:
    config = RunConfig(thermo={"population_scope": "all_approximate"})
    records = [
        _thermo("a", "mol1", "C2H6O", 6, 0.0),
        _thermo("b", "mol2", "C2H7N", 7, 1.0),
    ]

    ranked = compute_delta_g_and_populations(records, config)

    assert all(record.approximate_population for record in ranked)
    assert all(CROSS_PROTOMER_WARNING in record.warnings for record in ranked)


def test_population_correction_provider_can_shift_mixed_microstates() -> None:
    class MockCorrectionProvider:
        def get_microstate_correction(self, record, ph, solvent, temperature):
            assert ph == 7.0
            assert solvent == "water"
            assert temperature == pytest.approx(298.15)
            return {"a": 0.0, "b": -2.0}[record.id]

    config = RunConfig(thermo={"population_scope": "same_charge"})
    records = [
        _thermo("a", "mol1", "C2H6O", 6, 0.0, charge=0),
        _thermo("b", "mol2", "C2H7N", 7, 1.0, charge=0),
    ]

    ranked = compute_delta_g_and_populations(
        records,
        config,
        correction_provider=MockCorrectionProvider(),
    )

    assert ranked[0].parent_id == "b"
    assert ranked[0].relative_free_energy_kcal_mol == 0.0
    assert ranked[0].score_kcal_mol == -1.0
    assert not ranked[0].approximate_population
    assert CROSS_PROTOMER_WARNING not in ranked[0].warnings
    assert ranked[0].metadata["ranking"]["microstate_correction_kcal_mol"] == -2.0


def test_ranked_outputs_are_written(tmp_path: Path) -> None:
    config = RunConfig(thermo={"population_scope": "same_formula"})
    ranked = compute_delta_g_and_populations(
        [_thermo("a", "mol", "C2H6O", 6, 0.0)],
        config,
    )

    write_ranked_outputs(ranked, tmp_path)

    assert (tmp_path / "ranked_variants.csv").exists()
    assert "population_is_approximate" in (
        tmp_path / "ranked_variants.csv"
    ).read_text(encoding="utf-8").splitlines()[0]
    assert (tmp_path / "ranked_variants.sdf").exists()
    assert (tmp_path / "ranked_variants.json").exists()
    assert (tmp_path / "ranking_summary.md").exists()


def test_run_xtb_thermo_with_mocked_subprocess(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    conformer = _crest(tmp_path)
    config = RunConfig(input_path=tmp_path / "input.smi", output_dir=tmp_path / "run")
    monkeypatch.setattr(xtb_runner.shutil, "which", lambda executable: f"/mock/{executable}")

    def fake_run_command(command, cwd=None, **kwargs):
        assert cwd is not None
        workdir = Path(cwd)
        if "thermo" in command:
            (workdir / "xtb_thermo.out").write_text(
                Path("tests/data/xtb_thermo.log").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
        return subprocess.CompletedProcess(command, 0, stdout="xtb ok\n", stderr="")

    monkeypatch.setattr(xtb_runner, "run_command", fake_run_command)

    record = run_xtb_thermo(conformer, config)

    assert record.parent_id == conformer.id
    assert record.free_energy_kcal_mol == pytest.approx(hartree_to_kcal_mol(-39.98))
    assert record.entropy_cal_mol_k == 21.5
    assert (tmp_path / "run").glob("**/xtb_thermo.json")


def _thermo(
    record_id: str,
    molname: str,
    formula: str,
    proton_count: int,
    free_energy: float,
    charge: int = 0,
) -> ThermoRecord:
    return ThermoRecord(
        id=record_id,
        parent_id=f"crest_{record_id}",
        input_molecule_id=molname,
        molname=molname,
        canonical_smiles="CCO",
        isomeric_smiles="CCO",
        molecular_formula=formula,
        formal_charge=charge,
        explicit_proton_count=proton_count,
        source_software="test",
        source_python_function="test",
        free_energy_kcal_mol=free_energy,
    )


def _crest(tmp_path: Path) -> CrestConformerRecord:
    workdir = tmp_path / "crest_workdir"
    workdir.mkdir()
    xyz_path = workdir / "crest_conformer_0001.xyz"
    xyz_path.write_text(
        "3\nmock\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 0.0 1.0 0.0\n",
        encoding="utf-8",
    )
    return CrestConformerRecord(
        id="crest_a",
        parent_id="seed_a",
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles="O",
        isomeric_smiles="O",
        molecular_formula="H2O",
        formal_charge=0,
        explicit_proton_count=2,
        source_software="crest",
        metadata={"crest": {"workdir": str(workdir)}},
        crest_index=1,
        energy_kcal_mol=0.0,
    )
