from pathlib import Path

from rdkit import Chem
from typer.testing import CliRunner

from dsvr.chemistry import conformers_auto3d
from dsvr.chemistry.conformers_auto3d import generate_auto3d_seeds
from dsvr.cli import app
from dsvr.config import RunConfig
from dsvr.models import (
    StereoRecord,
    make_input_id,
    make_protomer_id,
    make_stereo_id,
    make_tautomer_id,
)
from dsvr.runners.auto3d_runner import Auto3DUnavailableError


def test_generate_auto3d_seeds_from_mock_output_preserves_parent(
    tmp_path: Path,
    monkeypatch,
) -> None:
    stereo = _stereo("ethanol", "CCO")
    config = RunConfig(
        input_path=tmp_path / "stereo.sdf",
        output_dir=tmp_path / "run",
        seeding={"method": "auto3d", "auto3d_k": 2, "auto3d_model": "AIMNet2"},
    )

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
    ) -> tuple[Path, list[str]]:
        assert input_path.exists()
        assert k == 2
        assert model == "AIMNet2"
        assert internal_tautomer_stereo_enum is False
        output_sdf = output_dir / "mock_auto3d.sdf"
        _write_auto3d_output(output_sdf, "CCO", stereo.id, energy="-12.34")
        return output_sdf, ["auto3d", "run", str(input_path)]

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", fake_run_auto3d)

    records = generate_auto3d_seeds([stereo], config)

    assert len(records) == 1
    assert records[0].parent_id == stereo.id
    assert records[0].metadata["auto3d"]["lineage_mode"] == "post_stereo_seed"
    assert records[0].energy_kcal_mol == -12.34
    assert records[0].forcefield_status == "auto3d_optimized"
    assert any("disabled to avoid double enumeration" in w for w in records[0].warnings)
    seed_dir = tmp_path / "run" / "seeding" / "auto3d"
    assert (seed_dir / "auto3d_input.sdf").exists()
    assert (seed_dir / "auto3d_seeds.sdf").exists()
    assert (seed_dir / "auto3d_seeds.csv").exists()


def test_auto3d_internal_enum_marks_less_controlled_lineage(
    tmp_path: Path,
    monkeypatch,
) -> None:
    stereo = _stereo("ethanol", "CCO")
    config = RunConfig(
        input_path=tmp_path / "stereo.sdf",
        output_dir=tmp_path / "run",
        seeding={
            "method": "auto3d",
            "auto3d_internal_tautomer_stereo_enum": True,
        },
    )

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
    ) -> tuple[Path, list[str]]:
        assert internal_tautomer_stereo_enum is True
        output_sdf = output_dir / "mock_auto3d_internal.sdf"
        _write_auto3d_output(output_sdf, "CCO", None, energy="-1.0")
        return output_sdf, ["auto3d", "run", str(input_path)]

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", fake_run_auto3d)

    records = generate_auto3d_seeds([stereo], config)

    assert len(records) == 1
    assert records[0].parent_id == stereo.id
    assert records[0].metadata["auto3d"]["lineage_mode"] == "auto3d_internal_enum"
    assert any("less controlled" in warning for warning in records[0].warnings)


def test_cli_seed_auto3d_reports_missing_auto3d(tmp_path: Path, monkeypatch) -> None:
    stereo = _stereo("ethanol", "CCO")
    stereo_sdf = tmp_path / "stereo.sdf"
    _write_stereo_sdf(stereo_sdf, stereo)

    def missing_auto3d(*args, **kwargs):
        raise Auto3DUnavailableError("Install Auto3D before running seed-auto3d")

    monkeypatch.setattr("dsvr.cli.generate_auto3d_seeds", missing_auto3d)

    result = CliRunner().invoke(
        app,
        ["seed-auto3d", str(stereo_sdf), "--out", str(tmp_path / "out"), "--k", "5"],
    )

    assert result.exit_code != 0
    assert "Install Auto3D" in result.output


def test_cli_seed_auto3d_success_with_mocked_runner(tmp_path: Path, monkeypatch) -> None:
    stereo = _stereo("ethanol", "CCO")
    stereo_sdf = tmp_path / "stereo.sdf"
    _write_stereo_sdf(stereo_sdf, stereo)

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
    ) -> tuple[Path, list[str]]:
        assert internal_tautomer_stereo_enum is False
        output_sdf = output_dir / "mock_auto3d.sdf"
        _write_auto3d_output(output_sdf, "CCO", stereo.id, energy="-12.34")
        return output_sdf, ["auto3d", "run", str(input_path)]

    monkeypatch.setattr(conformers_auto3d, "run_auto3d", fake_run_auto3d)

    result = CliRunner().invoke(
        app,
        ["seed-auto3d", str(stereo_sdf), "--out", str(tmp_path / "out"), "--k", "5"],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "seeding" / "auto3d" / "auto3d_report.json").exists()
    assert "internal tautomer/stereo enumeration disabled" in result.output


def _stereo(molname: str, smiles: str) -> StereoRecord:
    molecule = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    input_id = make_input_id(molname, canonical)
    protomer_id = make_protomer_id(input_id, 1, canonical, isomeric)
    tautomer_id = make_tautomer_id(protomer_id, 1, canonical, isomeric)
    stereo_id = make_stereo_id(tautomer_id, 1, canonical, isomeric)
    return StereoRecord(
        id=stereo_id,
        parent_id=tautomer_id,
        input_molecule_id=input_id,
        molname=molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula="",
        formal_charge=Chem.GetFormalCharge(molecule),
        explicit_proton_count=None,
        source_software="test",
        source_python_function="test",
        stereo_index=1,
        rdkit_mol=molecule,
    )


def _write_auto3d_output(
    path: Path,
    smiles: str,
    stereo_id: str | None,
    *,
    energy: str,
) -> None:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    mol.SetProp("_Name", stereo_id or "auto3d_internal_output")
    if stereo_id is not None:
        mol.SetProp("DSVR_STEREO_ID", stereo_id)
    mol.SetProp("E_tot", energy)
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def _write_stereo_sdf(path: Path, stereo: StereoRecord) -> None:
    mol = Chem.Mol(stereo.rdkit_mol)
    mol.SetProp("_Name", stereo.id)
    mol.SetProp("DSVR_STEREO_ID", stereo.id)
    mol.SetProp("DSVR_INPUT_ID", stereo.input_molecule_id)
    mol.SetProp("DSVR_PARENT_TAUTOMER_ID", stereo.parent_id or "")
    mol.SetProp("DSVR_MOLNAME", stereo.molname)
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
