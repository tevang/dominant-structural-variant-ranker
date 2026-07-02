from pathlib import Path

from dsvr.io.sdf import read_sdf


def test_reads_multimolecule_sdf_names(tmp_path: Path) -> None:
    path = tmp_path / "mols.sdf"
    path.write_text(
        "mol_a\n  test\n\n  0  0  0  0  0  0  0  0  0  0  0  0 V2000\nM  END\n$$$$\n"
        "mol_b\n  test\n\n  0  0  0  0  0  0  0  0  0  0  0  0 V2000\nM  END\n$$$$\n",
        encoding="utf-8",
    )

    molecules = read_sdf(path)

    assert [mol.name for mol in molecules] == ["mol_a", "mol_b"]
    assert all(mol.input_hash for mol in molecules)

