import json
from pathlib import Path

from dsvr.models import (
    MoleculeRecord,
    ProtomerRecord,
    StereoRecord,
    TautomerRecord,
    make_input_id,
    make_protomer_id,
    make_stereo_id,
    make_tautomer_id,
)
from dsvr.workflow.provenance import write_all_provenance_outputs


def test_same_input_gives_same_ids() -> None:
    first = make_input_id("Ligand 1", "CCO", {"source": "unit"})
    second = make_input_id("Ligand 1", "CCO", {"source": "unit"})

    assert first == second
    assert first.startswith("ligand_1_")


def test_different_stereoisomers_have_different_ids() -> None:
    input_id = make_input_id("lactic", "CC(O)C(=O)O")
    protomer_id = make_protomer_id(input_id, 1, "CC(O)C(=O)O", "CC(O)C(=O)O")
    tautomer_id = make_tautomer_id(protomer_id, 1, "CC(O)C(=O)O", "CC(O)C(=O)O")

    first = make_stereo_id(tautomer_id, 1, "CC(O)C(=O)O", "C[C@H](O)C(=O)O")
    second = make_stereo_id(tautomer_id, 2, "CC(O)C(=O)O", "C[C@@H](O)C(=O)O")

    assert first != second


def test_parent_child_lineage_can_be_reconstructed_from_jsonl(tmp_path: Path) -> None:
    input_id = make_input_id("ethanol", "CCO")
    protomer_id = make_protomer_id(input_id, 1, "CCO", "CCO")
    tautomer_id = make_tautomer_id(protomer_id, 1, "CCO", "CCO")
    stereo_id = make_stereo_id(tautomer_id, 1, "CCO", "CCO")
    records = [
        MoleculeRecord(
            id=input_id,
            parent_id=None,
            input_molecule_id=input_id,
            molname="ethanol",
            canonical_smiles="CCO",
            isomeric_smiles="CCO",
            molecular_formula="C2H6O",
            formal_charge=0,
            explicit_proton_count=6,
            source_software="rdkit",
            source_python_function="test",
        ),
        ProtomerRecord(
            id=protomer_id,
            parent_id=input_id,
            input_molecule_id=input_id,
            molname="ethanol",
            canonical_smiles="CCO",
            isomeric_smiles="CCO",
            molecular_formula="C2H6O",
            formal_charge=0,
            explicit_proton_count=6,
            source_software="molscrub",
            source_python_function="test",
            protomer_index=1,
        ),
        TautomerRecord(
            id=tautomer_id,
            parent_id=protomer_id,
            input_molecule_id=input_id,
            molname="ethanol",
            canonical_smiles="CCO",
            isomeric_smiles="CCO",
            molecular_formula="C2H6O",
            formal_charge=0,
            explicit_proton_count=6,
            source_software="rdkit",
            source_python_function="test",
            tautomer_index=1,
        ),
        StereoRecord(
            id=stereo_id,
            parent_id=tautomer_id,
            input_molecule_id=input_id,
            molname="ethanol",
            canonical_smiles="CCO",
            isomeric_smiles="CCO",
            molecular_formula="C2H6O",
            formal_charge=0,
            explicit_proton_count=6,
            source_software="rdkit",
            source_python_function="test",
            stereo_index=1,
        ),
    ]

    write_all_provenance_outputs(records, tmp_path)

    lines = (tmp_path / "enumeration_provenance.jsonl").read_text(encoding="utf-8").splitlines()
    by_id = {record["id"]: record for record in (json.loads(line) for line in lines)}
    assert by_id[protomer_id]["parent_id"] == input_id
    assert by_id[tautomer_id]["parent_id"] == protomer_id
    assert by_id[stereo_id]["parent_id"] == tautomer_id
    assert (tmp_path / "inputs.csv").exists()
    assert (tmp_path / "protomers.csv").exists()
    assert (tmp_path / "tautomers.csv").exists()
    assert (tmp_path / "stereoisomers.csv").exists()
