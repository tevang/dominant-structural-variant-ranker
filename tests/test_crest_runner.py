import os
import subprocess
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from typer.testing import CliRunner

from dsvr import cli
from dsvr.config import RunConfig
from dsvr.models import (
    SeedConformerRecord,
    make_input_id,
    make_protomer_id,
    make_seed_id,
    make_stereo_id,
    make_tautomer_id,
)
from dsvr.parsing.crest_outputs import parse_crest_outputs
from dsvr.runners import crest_runner
from dsvr.runners.crest_runner import crest_workdir, read_seed_sdf, run_crest_for_seed


def test_parse_mocked_crest_outputs(tmp_path: Path) -> None:
    _write_mock_crest_outputs(tmp_path)

    parsed = parse_crest_outputs(tmp_path)

    assert len(parsed.conformers) == 2
    assert parsed.energy_source == tmp_path / "crest.energies"
    assert parsed.conformers[0].relative_energy_kcal_mol == 0.0
    assert parsed.conformers[1].relative_energy_kcal_mol == pytest.approx(1.255018948)


def test_run_crest_for_seed_with_mocked_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seed = _seed("ethanol", "CCO")
    config = RunConfig(input_path=tmp_path / "seeds.sdf", output_dir=tmp_path / "run")
    monkeypatch.setattr(crest_runner.shutil, "which", lambda executable: f"/mock/{executable}")

    def fake_run_command(command, cwd=None, **kwargs):
        assert cwd is not None
        if "--help" in command or "-h" in command:
            return subprocess.CompletedProcess(
                command,
                0,
                stdout="--gfn2 --chrg --uhf --alpb --gbsa --ewin -T\n",
                stderr="",
            )
        _write_mock_crest_outputs(Path(cwd))
        return subprocess.CompletedProcess(command, 0, stdout="crest done\n", stderr="")

    monkeypatch.setattr(crest_runner, "run_command", fake_run_command)

    records = run_crest_for_seed(seed, config)

    assert len(records) == 2
    assert all(record.parent_id == seed.id for record in records)
    assert records[0].input_molecule_id == seed.input_molecule_id
    assert records[0].source_software == "crest"
    assert records[0].energy_kcal_mol is not None
    assert records[0].relative_energy_kcal_mol == 0.0
    assert records[0].metadata["lineage"]["seed_id"] == seed.id
    assert records[0].metadata["lineage"]["stereo_id"] == seed.parent_id
    assert records[0].metadata["lineage"]["tautomer_id"] in seed.parent_id
    assert records[0].metadata["lineage"]["protomer_id"] in seed.parent_id
    workdir = crest_workdir(seed, config)
    assert workdir == tmp_path / "run" / "crest" / seed.input_molecule_id / seed.parent_id / seed.id
    assert (workdir / "input.xyz").exists()
    assert (workdir / "crest_provenance.jsonl").exists()
    assert (workdir / "crest_conformers.csv").exists()


def test_run_crest_failure_is_captured_in_json_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seed = _seed("ethanol", "CCO")
    config = RunConfig(
        input_path=tmp_path / "seeds.sdf",
        output_dir=tmp_path / "run",
        crest={"command_template": "crest {input_xyz} --mock-fail"},
    )

    def fake_run_command(command, cwd=None, **kwargs):
        raise crest_runner.ExternalToolError(
            "CREST failed",
            metadata={"returncode": 2, "failure_markers": ["error"]},
        )

    monkeypatch.setattr(crest_runner, "run_command", fake_run_command)

    records = run_crest_for_seed(seed, config)

    assert len(records) == 1
    assert records[0].crest_index == 0
    assert any("CREST failed" in warning for warning in records[0].warnings)
    workdir = crest_workdir(seed, config)
    assert (workdir / "crest_failures.json").exists()
    assert (workdir / "crest_provenance.jsonl").exists()


def test_read_seed_sdf_preserves_seed_lineage(tmp_path: Path) -> None:
    seed = _seed("ethanol", "CCO")
    seed_sdf = tmp_path / "seeds.sdf"
    _write_seed_sdf(seed_sdf, seed)

    records = read_seed_sdf(seed_sdf)

    assert len(records) == 1
    assert records[0].id == seed.id
    assert records[0].parent_id == seed.parent_id
    assert records[0].input_molecule_id == seed.input_molecule_id
    assert records[0].rdkit_mol is not None


def test_cli_run_crest_with_mocked_runner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seed = _seed("ethanol", "CCO")
    seed_sdf = tmp_path / "seeds.sdf"
    _write_seed_sdf(seed_sdf, seed)

    def fake_run_crest_for_seed(seed_record, config):
        workdir = config.output_dir / "crest" / seed_record.input_molecule_id
        workdir.mkdir(parents=True, exist_ok=True)
        return []

    monkeypatch.setattr(cli, "run_crest_for_seed", fake_run_crest_for_seed)

    result = CliRunner().invoke(
        cli.app,
        [
            "run-crest",
            str(seed_sdf),
            "--out",
            str(tmp_path / "out"),
            "--solvent",
            "water",
            "--ph",
            "7.0",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "crest" / "crest_report.json").exists()


@pytest.mark.skipif(
    os.environ.get("DSVR_RUN_EXTERNAL") != "1",
    reason="Real CREST execution is skipped unless DSVR_RUN_EXTERNAL=1.",
)
def test_real_crest_execution_smoke(tmp_path: Path) -> None:
    seed = _seed("ethanol", "CCO")
    config = RunConfig(input_path=tmp_path / "seeds.sdf", output_dir=tmp_path / "run")
    records = run_crest_for_seed(seed, config)
    assert records


def _seed(molname: str, smiles: str) -> SeedConformerRecord:
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(molecule, randomSeed=7)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    input_id = make_input_id(molname, canonical)
    protomer_id = make_protomer_id(input_id, 1, canonical, isomeric)
    tautomer_id = make_tautomer_id(protomer_id, 1, canonical, isomeric)
    stereo_id = make_stereo_id(tautomer_id, 1, canonical, isomeric)
    seed_id = make_seed_id(stereo_id, 1, canonical, isomeric)
    return SeedConformerRecord(
        id=seed_id,
        parent_id=stereo_id,
        input_molecule_id=input_id,
        molname=molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula="C2H6O",
        formal_charge=Chem.GetFormalCharge(molecule),
        explicit_proton_count=6,
        source_software="test",
        source_python_function="test",
        conformer_index=1,
        rdkit_mol=molecule,
        rdkit_conformer_id=0,
        embedding_status="success",
    )


def _write_mock_crest_outputs(workdir: Path) -> None:
    (workdir / "crest_conformers.xyz").write_text(
        "\n".join(
            [
                "3",
                "conformer 1",
                "O 0.0 0.0 0.0",
                "H 0.0 0.0 1.0",
                "H 0.0 1.0 0.0",
                "3",
                "conformer 2",
                "O 0.0 0.0 0.1",
                "H 0.0 0.0 1.1",
                "H 0.0 1.1 0.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (workdir / "crest.energies").write_text(
        "1 -40.000000\n2 -39.998000\n",
        encoding="utf-8",
    )


def _write_seed_sdf(path: Path, seed: SeedConformerRecord) -> None:
    mol = Chem.Mol(seed.rdkit_mol)
    mol.SetProp("_Name", seed.id)
    mol.SetProp("DSVR_SEED_ID", seed.id)
    mol.SetProp("DSVR_INPUT_ID", seed.input_molecule_id)
    mol.SetProp("DSVR_PARENT_STEREO_ID", seed.parent_id or "")
    mol.SetProp("DSVR_MOLNAME", seed.molname)
    mol.SetProp("DSVR_CANONICAL_SMILES", seed.canonical_smiles or "")
    mol.SetProp("DSVR_ISOMERIC_SMILES", seed.isomeric_smiles or "")
    mol.SetProp("DSVR_FORMULA", seed.molecular_formula or "")
    mol.SetProp("DSVR_FORMAL_CHARGE", str(seed.formal_charge))
    mol.SetProp("DSVR_EXPLICIT_PROTON_COUNT", str(seed.explicit_proton_count))
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
