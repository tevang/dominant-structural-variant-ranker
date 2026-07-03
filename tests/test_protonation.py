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
    ) -> tuple[list[Chem.Mol], str, str]:
        assert ph_low == 7.0
        assert ph_high == 7.0
        return [
            Chem.MolFromSmiles("CCN"),
            Chem.MolFromSmiles("CC[NH3+]"),
            Chem.MolFromSmiles("CC[NH3+]"),
        ], "molscrub-test", "Scrub(ph_low=7.0, ph_high=7.0)"

    monkeypatch.setattr(protonation, "generate_molscrub_candidates", fake_molscrub_candidates)
    config = RunConfig(
        input_path=input_path,
        output_dir=tmp_path / "run",
        enumeration={"max_protomers_per_molecule": 32},
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
    assert "candidate generation only" in records[0].warnings[0]
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
    ) -> tuple[list[Chem.Mol], str, str]:
        return [
            Chem.MolFromSmiles("CCN"),
            Chem.MolFromSmiles("CC[NH3+]"),
        ], "molscrub-test", "Scrub(...)"

    monkeypatch.setattr(protonation, "generate_molscrub_candidates", fake_molscrub_candidates)
    config = RunConfig(
        input_path=input_path,
        output_dir=tmp_path / "run",
        enumeration={"max_protomers_per_molecule": 1},
    )

    records = generate_protomer_candidates(molecules[0], config)

    assert len(records) == 1
    assert any("truncated to 1" in warning for warning in records[0].warnings)


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
