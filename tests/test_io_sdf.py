from pathlib import Path

from rdkit import Chem

from dsvr.io.read_inputs import read_molecules
from dsvr.io.sdf import read_sdf


def test_reads_multimolecule_sdf_names_and_properties(tmp_path: Path) -> None:
    path = tmp_path / "mols.sdf"
    writer = Chem.SDWriter(str(path))
    mol_a = Chem.MolFromSmiles("CCO")
    mol_a.SetProp("_Name", "mol_a")
    mol_a.SetProp("PUBCHEM_CID", "702")
    mol_b = Chem.MolFromSmiles("c1ccccc1")
    mol_b.SetProp("_Name", "mol_b")
    mol_b.SetProp("PUBCHEM_CID", "241")
    writer.write(mol_a)
    writer.write(mol_b)
    writer.close()

    molecules, invalid = read_sdf(path)

    assert invalid == []
    assert [mol.molname for mol in molecules] == ["mol_a", "mol_b"]
    assert molecules[0].input_properties["PUBCHEM_CID"] == "702"
    assert molecules[1].input_properties["PUBCHEM_CID"] == "241"
    assert all(mol.rdkit_mol is not None for mol in molecules)


def test_sdf_missing_name_generates_name(tmp_path: Path) -> None:
    path = tmp_path / "mols.sd"
    writer = Chem.SDWriter(str(path))
    mol = Chem.MolFromSmiles("CCO")
    writer.write(mol)
    writer.close()

    molecules = read_molecules(path)

    assert len(molecules) == 1
    assert molecules[0].molname == "mol_000001"
