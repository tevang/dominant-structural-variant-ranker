from pathlib import Path

from rdkit import Chem
from typer.testing import CliRunner

from dsvr.chemistry.stereochemistry import enumerate_stereoisomers
from dsvr.cli import app
from dsvr.config import RunConfig
from dsvr.models import TautomerRecord, make_input_id, make_protomer_id, make_tautomer_id


def test_enumerates_one_unassigned_chiral_center(tmp_path: Path) -> None:
    tautomer = _tautomer("lactic", "CC(O)C(=O)O")
    config = RunConfig(input_path=tmp_path / "input.sdf", output_dir=tmp_path / "run")

    records = enumerate_stereoisomers(tautomer, config)

    assert len(records) == 2
    assert {record.isomeric_smiles for record in records} == {
        "C[C@H](O)C(=O)O",
        "C[C@@H](O)C(=O)O",
    }
    assert all(record.parent_id == tautomer.id for record in records)


def test_enumerates_e_z_double_bond_stereochemistry(tmp_path: Path) -> None:
    tautomer = _tautomer("butene", "CC=CC")
    config = RunConfig(input_path=tmp_path / "input.sdf", output_dir=tmp_path / "run")

    records = enumerate_stereoisomers(tautomer, config)

    assert len(records) == 2
    assert {record.isomeric_smiles for record in records} == {"C/C=C/C", "C/C=C\\C"}


def test_stereoisomer_max_isomers_cap(tmp_path: Path) -> None:
    tautomer = _tautomer("lactic", "CC(O)C(=O)O")
    config = RunConfig(
        input_path=tmp_path / "input.sdf",
        output_dir=tmp_path / "run",
        enumeration={"max_stereoisomers_per_tautomer": 1},
    )

    records = enumerate_stereoisomers(tautomer, config)

    assert len(records) == 1
    assert any("max_stereoisomers_per_tautomer" in warning for warning in records[0].warnings)


def test_preserves_assigned_stereochemistry_by_default(tmp_path: Path) -> None:
    tautomer = _tautomer("assigned_lactic", "C[C@H](O)C(=O)O")
    config = RunConfig(input_path=tmp_path / "input.sdf", output_dir=tmp_path / "run")

    records = enumerate_stereoisomers(tautomer, config)

    assert len(records) == 1
    assert records[0].isomeric_smiles == "C[C@H](O)C(=O)O"
    assert any(
        "Assigned stereochemistry was preserved" in warning for warning in records[0].warnings
    )


def test_cli_enumerate_stereo_from_tautomer_sdf(tmp_path: Path) -> None:
    tautomer = _tautomer("lactic", "CC(O)C(=O)O")
    tautomers_sdf = tmp_path / "tautomers.sdf"
    writer = Chem.SDWriter(str(tautomers_sdf))
    mol = Chem.Mol(tautomer.rdkit_mol)
    mol.SetProp("_Name", tautomer.id)
    mol.SetProp("DSVR_TAUTOMER_ID", tautomer.id)
    mol.SetProp("DSVR_INPUT_ID", tautomer.input_molecule_id)
    mol.SetProp("DSVR_PARENT_PROTOMER_ID", tautomer.parent_id or "")
    mol.SetProp("DSVR_MOLNAME", tautomer.molname)
    writer.write(mol)
    writer.close()
    outdir = tmp_path / "out"

    result = CliRunner().invoke(
        app,
        [
            "enumerate-stereo",
            str(tautomers_sdf),
            "--out",
            str(outdir),
            "--max-isomers",
            "4",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "tryEmbedding is heuristic" in result.output
    assert (outdir / "enumeration" / "stereoisomers" / "stereo_report.json").exists()


def _tautomer(molname: str, smiles: str) -> TautomerRecord:
    molecule = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    input_id = make_input_id(molname, canonical)
    protomer_id = make_protomer_id(input_id, 1, canonical, isomeric)
    tautomer_id = make_tautomer_id(protomer_id, 1, canonical, isomeric)
    return TautomerRecord(
        id=tautomer_id,
        parent_id=protomer_id,
        input_molecule_id=input_id,
        molname=molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula="",
        formal_charge=0,
        explicit_proton_count=None,
        source_software="test",
        source_python_function="test",
        tautomer_index=1,
        rdkit_mol=molecule,
    )



def test_stereo_timeout_fallback_keeps_input_state(tmp_path: Path, monkeypatch) -> None:
    import time

    tautomer = _tautomer("lactic", "CC(O)C(=O)O")

    def hanging_worker(*args, **kwargs):
        time.sleep(10)

    monkeypatch.setattr("dsvr.chemistry.stereochemistry._stereo_worker", hanging_worker)
    config = RunConfig(
        input_path=tmp_path / "input.sdf",
        output_dir=tmp_path / "run",
        stereoisomer_filtering={"timeout_seconds_per_tautomer": 1},
    )

    records = enumerate_stereoisomers(tautomer, config)

    assert len(records) == 1
    assert records[0].isomeric_smiles == "CC(O)C(=O)O"
    assert any("STEREO_TIMEOUT_FALLBACK" in warning for warning in records[0].warnings)


def test_assigned_stereo_is_preserved_unless_assigned_centers_are_enabled(tmp_path: Path) -> None:
    tautomer = _tautomer("bromochlorofluoromethane", "F[C@H](Cl)Br")
    default_config = RunConfig(input_path=tmp_path / "input.sdf", output_dir=tmp_path / "default")
    default_records = enumerate_stereoisomers(tautomer, default_config)

    enumerate_assigned_config = RunConfig(
        input_path=tmp_path / "input.sdf",
        output_dir=tmp_path / "assigned",
        stereoisomer_filtering={"only_unassigned": False},
        enumeration={"stereo_only_unassigned": False},
    )
    assigned_records = enumerate_stereoisomers(tautomer, enumerate_assigned_config)

    assert len(default_records) == 1
    assert default_records[0].isomeric_smiles == "F[C@H](Cl)Br"
    assert {record.isomeric_smiles for record in assigned_records} == {
        "F[C@H](Cl)Br",
        "F[C@@H](Cl)Br",
    }
