from pathlib import Path

import pytest
from rdkit import Chem
from typer.testing import CliRunner

from dsvr import cli
from dsvr.chemistry import protonation
from dsvr.chemistry.protonation import generate_protomer_candidates
from dsvr.config import RunConfig
from dsvr.io.smiles import read_smiles
from dsvr.runners.molscrub_runner import MolscrubUnavailableError


def test_generate_protomer_candidates_with_mocked_molscrub(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCN ethylamine\n", encoding="utf-8")
    molecules, invalid = read_smiles(input_path)
    assert invalid == []

    def fake_molscrub_candidates(
        molecule: Chem.Mol,
        *,
        ph_low: float,
        ph_high: float,
        skip_gen3d: bool = True,
        timeout_seconds: int = 60,
    ) -> tuple[list[Chem.Mol], str, str]:
        assert ph_low == 7.0
        assert ph_high == 7.0
        assert skip_gen3d is True
        assert timeout_seconds == 60
        return [
            Chem.MolFromSmiles("CCN"),
            Chem.MolFromSmiles("CC[NH3+]"),
            Chem.MolFromSmiles("CC[NH3+]"),
        ], "molscrub-test", "Scrub(ph_low=7.0, ph_high=7.0)"

    monkeypatch.setattr(protonation, "generate_molscrub_candidates", fake_molscrub_candidates)
    config = RunConfig(
        input_path=input_path,
        output_dir=tmp_path / "run",
        protonation={"tool": "molscrub", "max_protomers_per_molecule": 32},
    )

    records = generate_protomer_candidates(molecules[0], config)

    assert len(records) == 2
    assert records[0].parent_id == molecules[0].input_id
    assert records[0].input_molecule_id == molecules[0].input_id
    assert records[0].molecular_formula == "C2H7N"
    assert records[0].formal_charge == 0
    assert records[0].explicit_proton_count == 7
    assert records[1].formal_charge == 1
    assert records[1].rdkit_mol is not None
    assert "candidate generation/filtering only" in records[0].warnings[0]
    protomer_dir = tmp_path / "run" / "enumeration" / "protomers"
    assert (protomer_dir / f"{molecules[0].input_id}_protomers.sdf").exists()
    assert (protomer_dir / f"{molecules[0].input_id}_protomers.csv").exists()


def test_generate_protomer_candidates_caps_and_warns(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCN ethylamine\n", encoding="utf-8")
    molecules, _ = read_smiles(input_path)

    def fake_molscrub_candidates(
        molecule: Chem.Mol,
        *,
        ph_low: float,
        ph_high: float,
        skip_gen3d: bool = True,
        timeout_seconds: int = 60,
    ) -> tuple[list[Chem.Mol], str, str]:
        return [
            Chem.MolFromSmiles("CCN"),
            Chem.MolFromSmiles("CC[NH3+]"),
        ], "molscrub-test", "Scrub(...)"

    monkeypatch.setattr(protonation, "generate_molscrub_candidates", fake_molscrub_candidates)
    config = RunConfig(
        input_path=input_path,
        output_dir=tmp_path / "run",
        protonation={"tool": "molscrub", "max_protomers_per_molecule": 1},
    )

    records = generate_protomer_candidates(molecules[0], config)

    assert len(records) == 1
    protomer_dir = tmp_path / "run" / "enumeration" / "protomers"
    rejected = protomer_dir / "protomers_rejected.csv"
    assert rejected.exists()
    assert "beyond_max_protomers_per_molecule" in rejected.read_text(encoding="utf-8") or "lower_scoring_same_charge_state" in rejected.read_text(encoding="utf-8")


def test_cli_enumerate_protomers_missing_molscrub_is_actionable(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCN ethylamine\n", encoding="utf-8")

    def missing_molscrub(*args: object, **kwargs: object) -> list[object]:
        raise MolscrubUnavailableError("install molscrub with pip install git+https://github.com/forlilab/molscrub.git")

    monkeypatch.setattr(cli, "generate_protomer_candidates", missing_molscrub)
    result = CliRunner().invoke(
        cli.app,
        [
            "enumerate-protomers",
            str(input_path),
            "--ph",
            "7.0",
            "--solvent",
            "water",
            "--out",
            str(tmp_path / "out"),
        ],
    )

    assert result.exit_code != 0
    assert "install molscrub" in result.output


def test_cli_enumerate_protomers_with_mocked_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCN ethylamine\n", encoding="utf-8")
    outdir = tmp_path / "out"

    def fake_generate(molecule: object, config: RunConfig) -> list[object]:
        protomer_dir = config.output_dir / "enumeration" / "protomers"
        protomer_dir.mkdir(parents=True, exist_ok=True)
        return []

    monkeypatch.setattr(cli, "generate_protomer_candidates", fake_generate)
    result = CliRunner().invoke(
        cli.app,
        [
            "enumerate-protomers",
            str(input_path),
            "--ph",
            "7.0",
            "--solvent",
            "water",
            "--out",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    report = outdir / "enumeration" / "protomers" / "protomer_report.json"
    assert report.exists()
    assert "candidate generation only" in result.output


def test_many_molscrub_states_only_selected_protomers_pass_downstream(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCN ethylamine\n", encoding="utf-8")
    molecules, _ = read_smiles(input_path)

    def fake_molscrub_candidates(
        molecule: Chem.Mol,
        *,
        ph_low: float,
        ph_high: float,
        skip_gen3d: bool = True,
        timeout_seconds: int = 60,
    ) -> tuple[list[Chem.Mol], str, str]:
        return [
            Chem.MolFromSmiles("CCN"),
            Chem.MolFromSmiles("CC[NH3+]"),
            Chem.MolFromSmiles("C[NH2+]C"),
            Chem.MolFromSmiles("C[NH+]C"),
            Chem.MolFromSmiles("C[NH3+]C"),
            Chem.MolFromSmiles("CC[NH3+]"),
        ], "molscrub-test", "Scrub(skip_gen3d=True)"

    monkeypatch.setattr(protonation, "generate_molscrub_candidates", fake_molscrub_candidates)
    config = RunConfig(
        input_path=input_path,
        output_dir=tmp_path / "run",
        protonation={"tool": "molscrub", "max_protomers_per_molecule": 2},
    )

    records = generate_protomer_candidates(molecules[0], config)

    assert len(records) == 2
    assert all(record.metadata["plausibility"]["selected"] for record in records)
    protomer_dir = tmp_path / "run" / "enumeration" / "protomers"
    assert (protomer_dir / "protomers_all.csv").exists()
    assert (protomer_dir / "protomers_selected.csv").exists()
    assert (protomer_dir / "protomers_rejected.csv").exists()
    assert "score_is_population_estimate" in (protomer_dir / "protomers_selected.csv").read_text(encoding="utf-8")


def test_no_valid_molscrub_state_retains_original_with_warning(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCN ethylamine\n", encoding="utf-8")
    molecules, _ = read_smiles(input_path)

    def fake_molscrub_candidates(
        molecule: Chem.Mol,
        *,
        ph_low: float,
        ph_high: float,
        skip_gen3d: bool = True,
        timeout_seconds: int = 60,
    ) -> tuple[list[Chem.Mol], str, str]:
        return [], "molscrub-test", "Scrub(...)"

    monkeypatch.setattr(protonation, "generate_molscrub_candidates", fake_molscrub_candidates)
    config = RunConfig(
        input_path=input_path,
        output_dir=tmp_path / "run",
        protonation={"tool": "molscrub"},
    )

    records = generate_protomer_candidates(molecules[0], config)

    assert len(records) == 1
    assert any("retained input" in warning for warning in records[0].warnings)


def _unipka_molecule(tmp_path: Path, smiles: str):
    input_path = tmp_path / "mols.smi"
    input_path.write_text(f"{smiles} testmol\n", encoding="utf-8")
    molecules, invalid = read_smiles(input_path)
    assert invalid == []
    return molecules[0]


def _unipka_result(input_id: str, forms, envelope=None, failed=False, warning=None):
    from dsvr.runners.unipka_runner import UnipkaEnvelope, UnipkaForm, UnipkaMoleculeResult

    return UnipkaMoleculeResult(
        input_id=input_id,
        forms=[UnipkaForm(form_smiles=smi, occupancy=occ) for smi, occ in forms],
        envelope=UnipkaEnvelope(microstates=envelope or {}),
        failed=failed,
        failed_warning=warning,
        source_command="apptainer run container.sif protonate ...",
    )


def test_unipka_selection_keeps_dominant_form_above_threshold(tmp_path: Path) -> None:
    from dsvr.chemistry.protonation import generate_unipka_protomer_candidates

    molecule = _unipka_molecule(tmp_path, "CC(=O)O")
    result = _unipka_result(
        molecule.input_id,
        [("CC(=O)[O-]", 0.97), ("CC(=O)O", 0.03)],
        envelope={"CC(=O)[O-]": -6.79, "CC(=O)O": -5.25},
    )
    config = RunConfig(
        input_path=tmp_path / "mols.smi",
        output_dir=tmp_path / "run",
        protonation={"tool": "unipka", "max_protomers_per_molecule": 4},
    )

    records = generate_unipka_protomer_candidates(molecule, config, result)

    charges = [record.formal_charge for record in records]
    assert charges.count(-1) == 1
    assert len(records) == 2
    # keep_input_state has selection precedence over raw occupancy ordering,
    # so the neutral input form comes first; occupancy metadata matches each form
    by_smiles = {record.canonical_smiles: record for record in records}
    assert by_smiles["CC(=O)[O-]"].metadata["unipka_occupancy"] == 0.97
    assert by_smiles["CC(=O)[O-]"].metadata["unipka_dg"] == -6.79
    assert by_smiles["CC(=O)O"].metadata["unipka_occupancy"] == 0.03
    assert by_smiles["CC(=O)O"].metadata["plausibility"]["selection_reason"] == "unipka_occupancy_ranked"
    assert by_smiles["CC(=O)[O-]"].metadata["plausibility"]["selection_reason"] == "unipka_occupancy_ranked"


def test_unipka_cap_trims_by_occupancy(tmp_path: Path) -> None:
    from dsvr.chemistry.protonation import generate_unipka_protomer_candidates

    molecule = _unipka_molecule(tmp_path, "c1c[nH]cn1")
    result = _unipka_result(
        molecule.input_id,
        [
            ("c1c[nH]cn1", 0.60),
            ("c1c[nH+]c[nH]1", 0.35),
            ("c1c[nH]c[nH+]1", 0.05),
        ],
        envelope={},
    )
    config = RunConfig(
        input_path=tmp_path / "mols.smi",
        output_dir=tmp_path / "run",
        protonation={"tool": "unipka", "max_protomers_per_molecule": 2, "unipka": {"max_forms": 3}},
    )

    records = generate_unipka_protomer_candidates(molecule, config, result)

    assert len(records) == 2
    assert [r.metadata["unipka_occupancy"] for r in records] == [0.60, 0.35]


def test_unipka_failed_molecule_retains_input_state_with_summary_null(tmp_path: Path) -> None:
    from dsvr.chemistry.protonation import generate_unipka_protomer_candidates

    molecule = _unipka_molecule(tmp_path, "C1CCCCC1")
    result = _unipka_result(
        molecule.input_id,
        [],
        failed=True,
        warning="Uni-Pka returned no protonation form; retained input state",
    )
    config = RunConfig(
        input_path=tmp_path / "mols.smi",
        output_dir=tmp_path / "run",
        protonation={"tool": "unipka"},
    )

    records = generate_unipka_protomer_candidates(molecule, config, result)

    assert len(records) == 1
    assert any("retained input state" in warning for warning in records[0].warnings)
    summary = records[0].metadata["unipka_summary"]
    assert summary["microstate_count"] == 0
    assert summary["occupancy_entropy"] is None
    assert records[0].metadata["unipka_occupancy"] is None
    assert records[0].source_software == "unipka"


def test_unipka_summary_computed_from_envelope(tmp_path: Path) -> None:
    from dsvr.chemistry.protonation import generate_unipka_protomer_candidates

    molecule = _unipka_molecule(tmp_path, "c1c[nH]cn1")
    result = _unipka_result(
        molecule.input_id,
        [("c1c[nH]cn1", 0.6254), ("c1c[nH+]c[nH]1", 0.3746)],
        envelope={"c1c[nH]cn1": -5.24537, "c1c[nH+]c[nH]1": -6.7938},
    )
    config = RunConfig(
        input_path=tmp_path / "mols.smi",
        output_dir=tmp_path / "run",
        chemistry={"ph": 7.4},
        protonation={"tool": "unipka"},
    )

    records = generate_unipka_protomer_candidates(molecule, config, result)

    summary = records[0].metadata["unipka_summary"]
    assert summary["microstate_count"] == 2
    assert summary["charge_population"] == {"0": pytest.approx(0.6254, abs=1e-3), "1": pytest.approx(0.3746, abs=1e-3)}
    assert summary["top_two_occupancy_gap"] == pytest.approx(0.6254 - 0.3746, abs=1e-3)
    # imidazole pKa ~7.18, working pH 7.4 → small nearest distance
    assert summary["pka_nearest_transition"] == pytest.approx(7.18, abs=0.2)
    # neutral↔+1 only: never negative → pI null
    assert summary["isoelectric_point"] is None


def test_unipka_provenance_records_tool_and_command(tmp_path: Path) -> None:
    from dsvr.chemistry.protonation import generate_unipka_protomer_candidates

    molecule = _unipka_molecule(tmp_path, "CC(=O)O")
    result = _unipka_result(
        molecule.input_id,
        [("CC(=O)[O-]", 0.97)],
        envelope={"CC(=O)[O-]": -6.79},
    )
    config = RunConfig(
        input_path=tmp_path / "mols.smi",
        output_dir=tmp_path / "run",
        protonation={"tool": "unipka"},
    )

    records = generate_unipka_protomer_candidates(molecule, config, result)

    assert records[0].source_software == "unipka"
    assert "protonate" in (records[0].source_command or "")
    assert records[0].source_python_function and "unipka" in records[0].source_python_function
    assert records[0].metadata["tool"] == "unipka"
    assert records[0].metadata["target_ph"] == 7.0
