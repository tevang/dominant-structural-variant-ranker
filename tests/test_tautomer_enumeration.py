from pathlib import Path

from rdkit import Chem
from typer.testing import CliRunner

from dsvr.chemistry.tautomers import enumerate_tautomers
from dsvr.cli import app
from dsvr.config import RunConfig
from dsvr.models import ProtomerRecord, make_input_id, make_protomer_id


def test_enumerate_tautomers_known_example(tmp_path: Path) -> None:
    protomer = _protomer("acetylacetone", "CC(=O)CC(C)=O")
    config = RunConfig(
        input_path=tmp_path / "input.sdf",
        output_dir=tmp_path / "run",
        enumeration={"max_tautomers_per_protomer": 16},
    )

    records = enumerate_tautomers(protomer, config)

    assert len(records) >= 2
    assert all(record.parent_id == protomer.id for record in records)
    assert all(record.input_molecule_id == protomer.input_molecule_id for record in records)
    assert all(record.rdkit_mol is not None for record in records)
    tautomer_dir = tmp_path / "run" / "enumeration" / "tautomers"
    assert (tautomer_dir / f"{protomer.id}_tautomers.sdf").exists()
    assert (tautomer_dir / f"{protomer.id}_tautomers.csv").exists()


def test_tautomer_records_do_not_claim_stability_ranking(tmp_path: Path) -> None:
    protomer = _protomer("acetylacetone", "CC(=O)CC(C)=O")
    config = RunConfig(input_path=tmp_path / "input.sdf", output_dir=tmp_path / "run")

    records = enumerate_tautomers(protomer, config)

    assert records
    assert all(record.metadata["candidate_generation_only"] is True for record in records)
    assert all(record.metadata["not_stability_ranking"] is True for record in records)
    assert any("no tautomer stability ranking is implied" in records[0].warnings[0] for _ in [0])


def test_tautomer_enumeration_cap_warns(tmp_path: Path) -> None:
    protomer = _protomer("acetylacetone", "CC(=O)CC(C)=O")
    config = RunConfig(
        input_path=tmp_path / "input.sdf",
        output_dir=tmp_path / "run",
        enumeration={"max_tautomers_per_protomer": 1},
    )

    records = enumerate_tautomers(protomer, config)

    assert len(records) == 1
    assert any("max_tautomers_per_protomer" in warning for warning in records[0].warnings)


def test_cli_enumerate_tautomers_from_protomer_sdf(tmp_path: Path) -> None:
    protomer = _protomer("acetylacetone", "CC(=O)CC(C)=O")
    protomers_sdf = tmp_path / "protomers.sdf"
    writer = Chem.SDWriter(str(protomers_sdf))
    mol = Chem.Mol(protomer.rdkit_mol)
    mol.SetProp("_Name", protomer.id)
    mol.SetProp("DSVR_PROTOMER_ID", protomer.id)
    mol.SetProp("DSVR_INPUT_ID", protomer.input_molecule_id)
    mol.SetProp("DSVR_PARENT_ID", protomer.parent_id or "")
    mol.SetProp("DSVR_MOLNAME", protomer.molname)
    writer.write(mol)
    writer.close()
    outdir = tmp_path / "out"

    result = CliRunner().invoke(
        app,
        [
            "enumerate-tautomers",
            str(protomers_sdf),
            "--out",
            str(outdir),
            "--max-tautomers",
            "4",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "not tautomer stability ranking" in result.output
    assert (outdir / "enumeration" / "tautomers" / "tautomer_report.json").exists()


def _protomer(molname: str, smiles: str) -> ProtomerRecord:
    molecule = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    input_id = make_input_id(molname, canonical)
    protomer_id = make_protomer_id(input_id, 1, canonical, isomeric)
    return ProtomerRecord(
        id=protomer_id,
        parent_id=input_id,
        input_molecule_id=input_id,
        molname=molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula="C5H8O2",
        formal_charge=0,
        explicit_proton_count=8,
        source_software="test",
        source_python_function="test",
        protomer_index=1,
        rdkit_mol=molecule,
    )
