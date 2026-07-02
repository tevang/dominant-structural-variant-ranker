from pathlib import Path

from dsvr.io.smiles import read_smiles


def test_reads_two_column_smiles(tmp_path: Path) -> None:
    path = tmp_path / "mols.smi"
    path.write_text("CCO ethanol\nc1ccccc1 benzene\n", encoding="utf-8")

    molecules = read_smiles(path)

    assert [mol.smiles for mol in molecules] == ["CCO", "c1ccccc1"]
    assert [mol.name for mol in molecules] == ["ethanol", "benzene"]

