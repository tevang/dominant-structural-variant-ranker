from pathlib import Path

from dsvr.io.read_inputs import read_molecules, validate_input_file
from dsvr.io.smiles import read_smiles


def test_reads_two_column_smiles_with_names(tmp_path: Path) -> None:
    path = tmp_path / "mols.smi"
    path.write_text("CCO ethanol\nc1ccccc1 benzene\n", encoding="utf-8")

    molecules, invalid = read_smiles(path)

    assert invalid == []
    assert [mol.original_smiles for mol in molecules] == ["CCO", "c1ccccc1"]
    assert [mol.molname for mol in molecules] == ["ethanol", "benzene"]
    assert all(mol.canonical_smiles for mol in molecules)
    assert all(mol.rdkit_mol is not None for mol in molecules)


def test_reads_smiles_without_names_generates_zero_padded_names(tmp_path: Path) -> None:
    path = tmp_path / "mols.smi"
    path.write_text("CCO\nCCN\n", encoding="utf-8")

    molecules, invalid = read_smiles(path)

    assert invalid == []
    assert [mol.molname for mol in molecules] == ["mol_000001", "mol_000002"]
    assert [mol.input_id for mol in molecules] == ["mol_000001", "mol_000002"]


def test_reads_tabular_smiles_header(tmp_path: Path) -> None:
    path = tmp_path / "mols.smi"
    path.write_text("SMILES molname\nCCO ethanol\nCCN ethylamine\n", encoding="utf-8")

    molecules, invalid = read_smiles(path)

    assert invalid == []
    assert [mol.molname for mol in molecules] == ["ethanol", "ethylamine"]


def test_invalid_smiles_are_written_to_csv(tmp_path: Path) -> None:
    path = tmp_path / "mols.smi"
    invalid_path = tmp_path / "invalid_inputs.csv"
    path.write_text("CCO ethanol\nnot_a_smiles bad\n", encoding="utf-8")

    molecules, invalid = validate_input_file(path, invalid_output_path=invalid_path)

    assert len(molecules) == 1
    assert len(invalid) == 1
    assert invalid_path.exists()
    assert "RDKit failed to parse SMILES" in invalid_path.read_text(encoding="utf-8")


def test_read_molecules_supports_txt_suffix(tmp_path: Path) -> None:
    path = tmp_path / "mols.txt"
    path.write_text("CCO ethanol\n", encoding="utf-8")

    molecules = read_molecules(path)

    assert len(molecules) == 1
    assert molecules[0].molname == "ethanol"


def test_parse_all_supplied_example_smiles() -> None:
    molecules, invalid = read_smiles(Path("examples/test_molecules.smi"))

    assert invalid == []
    assert len(molecules) == 8
    assert [molecule.molname for molecule in molecules] == [
        "4862293",
        "1544787",
        "133506781",
        "38898616",
        "65986444",
        "170222839",
        "68880434",
        "4838114",
    ]
    assert all(molecule.rdkit_mol is not None for molecule in molecules)
