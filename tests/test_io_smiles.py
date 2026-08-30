from pathlib import Path

from dsvr.io.read_inputs import INVALID_INPUT_COLUMNS, read_molecules, validate_input_file
from dsvr.io.sdf import read_sdf
from dsvr.io.smiles import read_smiles

INVALID_RECORD_KEYS = set(INVALID_INPUT_COLUMNS)


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


def test_read_molecules_supports_csv_suffix(tmp_path: Path) -> None:
    path = tmp_path / "mols.csv"
    path.write_text("smiles,name\nCCO,ethanol\n", encoding="utf-8")

    molecules = read_molecules(path)

    assert len(molecules) == 1
    assert molecules[0].molname == "ethanol"
    assert molecules[0].original_smiles == "CCO"


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


def test_reads_identifier_first_csv_with_recognized_headers(tmp_path: Path) -> None:
    path = tmp_path / "chembl.csv"
    path.write_text(
        "chembl_id,canonical_smiles\nCHEMBL1993996,Cc1ccccc1\nCHEMBL25,CCO\n",
        encoding="utf-8",
    )

    molecules, invalid = read_smiles(path)

    assert invalid == []
    assert [mol.original_smiles for mol in molecules] == ["Cc1ccccc1", "CCO"]
    assert [mol.molname for mol in molecules] == ["CHEMBL1993996", "CHEMBL25"]
    assert [mol.input_id for mol in molecules] == ["mol_000001", "mol_000002"]


def test_reads_name_smiles_header_order_with_semicolons_and_tabs(tmp_path: Path) -> None:
    semicolon = tmp_path / "mols_semicolon.csv"
    semicolon.write_text("name;smiles\nethanol;CCO\nbenzene;c1ccccc1\n", encoding="utf-8")
    tabbed = tmp_path / "mols_tab.smi"
    tabbed.write_text("name\tsmiles\nethanol\tCCO\n", encoding="utf-8")

    semi_molecules, semi_invalid = read_smiles(semicolon)
    tab_molecules, tab_invalid = read_smiles(tabbed)

    assert semi_invalid == [] and tab_invalid == []
    assert [mol.original_smiles for mol in semi_molecules] == ["CCO", "c1ccccc1"]
    assert [mol.molname for mol in semi_molecules] == ["ethanol", "benzene"]
    assert tab_molecules[0].original_smiles == "CCO"
    assert tab_molecules[0].molname == "ethanol"


def test_probe_resolves_unrecognized_header(tmp_path: Path) -> None:
    path = tmp_path / "export.csv"
    path.write_text(
        "Molecule ChEMBL ID,Structure\n"
        "CHEMBL1993996,Cc1ccccc1\n"
        "CHEMBL25,CCO\n"
        "CHEMBL521,CC(=O)O\n",
        encoding="utf-8",
    )

    molecules, invalid = read_smiles(path)

    assert invalid == []
    assert [mol.original_smiles for mol in molecules] == ["Cc1ccccc1", "CCO", "CC(=O)O"]
    assert [mol.molname for mol in molecules] == ["CHEMBL1993996", "CHEMBL25", "CHEMBL521"]


def test_probe_prefers_earliest_column_on_parse_ties(tmp_path: Path) -> None:
    path = tmp_path / "ambiguous.smi"
    path.write_text("CCO CCN\nCCN CCO\n", encoding="utf-8")

    molecules, invalid = read_smiles(path)

    assert invalid == []
    assert [mol.original_smiles for mol in molecules] == ["CCO", "CCN"]


def test_no_smiles_column_marks_every_record_invalid(tmp_path: Path) -> None:
    path = tmp_path / "ids_only.csv"
    path.write_text("foo,bar\nbaz,qux\nquux,corge\n", encoding="utf-8")

    molecules, invalid = read_smiles(path)

    assert molecules == []
    assert len(invalid) == 3
    assert all(record["error"] == "no SMILES column found" for record in invalid)
    assert all(record["smiles"] == "" for record in invalid)
    assert all(set(record) == INVALID_RECORD_KEYS for record in invalid)


def test_data_named_like_a_header_name_is_not_treated_as_header(tmp_path: Path) -> None:
    path = tmp_path / "named_id.smi"
    path.write_text("CCO id\nCCN id2\n", encoding="utf-8")

    molecules, invalid = read_smiles(path)

    assert invalid == []
    assert [mol.molname for mol in molecules] == ["id", "id2"]


def test_invalid_record_carries_fixed_fields(tmp_path: Path) -> None:
    path = tmp_path / "mols.smi"
    path.write_text("CCO ethanol\nnot_a_smiles bad\n", encoding="utf-8")

    _, invalid = read_smiles(path)

    assert len(invalid) == 1
    record = invalid[0]
    assert set(record) == INVALID_RECORD_KEYS
    assert record["input_id"] == "mol_000002"
    assert record["source_format"] == "smiles"
    assert record["line_number"] == "2"
    assert record["name"] == "bad"
    assert record["smiles"] == "not_a_smiles"
    assert record["raw_record"] == "not_a_smiles bad"
    assert record["error"] == "RDKit failed to parse SMILES"


def test_invalid_inputs_csv_replaced_on_every_run(tmp_path: Path) -> None:
    valid = tmp_path / "mols.smi"
    valid.write_text("CCO ethanol\n", encoding="utf-8")
    invalid_path = tmp_path / "invalid_inputs.csv"
    invalid_path.write_text("record_index,raw_record,error\n0,stale,boom\n", encoding="utf-8")

    molecules, invalid = validate_input_file(valid, invalid_output_path=invalid_path)

    assert len(molecules) == 1 and invalid == []
    content = invalid_path.read_text(encoding="utf-8")
    assert "stale" not in content
    assert content == ",".join(INVALID_INPUT_COLUMNS) + "\n"

    rejecting = tmp_path / "rejecting.smi"
    rejecting.write_text("not_a_smiles bad\n", encoding="utf-8")
    _, invalid = validate_input_file(rejecting, invalid_output_path=invalid_path)

    lines = invalid_path.read_text(encoding="utf-8").splitlines()
    assert invalid and len(invalid) == 1
    assert lines[0] == ",".join(INVALID_INPUT_COLUMNS)
    assert len(lines) == 2
    assert "not_a_smiles" in lines[1]


def test_sdf_invalid_records_share_the_fixed_field_set(tmp_path: Path) -> None:
    smiles_path = tmp_path / "mols.smi"
    smiles_path.write_text("not_a_smiles bad\n", encoding="utf-8")
    sdf_path = tmp_path / "broken.sdf"
    sdf_path.write_text("broken title\nnot a real sdf record\n$$$$\n", encoding="utf-8")

    _, smiles_invalid = read_smiles(smiles_path)
    _, sdf_invalid = read_sdf(sdf_path)

    assert smiles_invalid and sdf_invalid
    assert set(sdf_invalid[0]) == set(smiles_invalid[0]) == INVALID_RECORD_KEYS
    assert sdf_invalid[0]["source_format"] == "sdf"
    assert sdf_invalid[0]["name"] == "broken title"
