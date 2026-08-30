"""Tests for bounded stereo over-enumeration (fill-to-cap) and the
cross-tautomer exact-duplicate stereoisomer guard."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
from rdkit import Chem

from dsvr.chemistry import stereo_auto3d_filter
from dsvr.chemistry.stereo_auto3d_filter import filter_stereoisomers_with_auto3d
from dsvr.chemistry.stereochemistry import (
    STEREO_ENUMERATION_HARD_CEILING,
    _enumeration_ceiling,
    enumerate_stereoisomers,
)
from dsvr.config import RunConfig
from dsvr.models import (
    StereoRecord,
    TautomerRecord,
    make_input_id,
    make_protomer_id,
    make_tautomer_id,
)


def _tautomer(smiles: str = "CC(O)C(=O)O", molname: str = "mol") -> TautomerRecord:
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


def test_bounded_over_enumeration_fills_duplicate_rich_enumeration_to_cap(
    tmp_path: Path,
    monkeypatch,
) -> None:
    tautomer = _tautomer()
    smiles_by_label = {
        "A": "C[C@H](O)C(=O)O",
        "B": "C[C@@H](O)C(=O)O",
        "C": "C[C@H](O)C(=O)N",
        "D": "C[C@@H](O)C(=O)N",
    }
    duplicate_rich = [
        smiles_by_label[label] for label in ("A", "A", "B", "A", "C", "B", "D")
    ]
    captured: dict[str, int] = {}

    def fake_enumerate(molecule, *, max_isomers, **kwargs):
        captured["max_isomers"] = max_isomers
        return [Chem.MolFromSmiles(smiles) for smiles in duplicate_rich[:max_isomers]]

    monkeypatch.setattr(
        "dsvr.chemistry.stereochemistry._enumerate_with_timeout", fake_enumerate
    )
    config = RunConfig(
        input_path=tmp_path / "input.sdf",
        output_dir=tmp_path / "run",
        enumeration={"max_stereoisomers_per_tautomer": 3},
    )

    records = enumerate_stereoisomers(tautomer, config)

    # Bounded ceiling: cap * enumeration_ceiling_multiplier (3 * 4 = 12).
    assert captured["max_isomers"] == 12
    selected = [record for record in records if record.metadata.get("selected", True)]
    pool = [record for record in records if not record.metadata.get("selected", True)]
    # All 4 unique isomers exist within the raw enumeration, so the first
    # cap=3 unique are selected and the 4th stays in the refill pool.
    assert len(selected) == 3
    assert len({record.isomeric_smiles for record in selected}) == 3
    assert len(pool) == 1
    assert "beyond_max_stereoisomers_per_tautomer" in str(pool[0].metadata)
    # Pool records are not part of the per-tautomer SDF outputs.
    supplier = Chem.SDMolSupplier(
        str(tmp_path / "run" / "enumeration" / "stereoisomers" / f"{tautomer.id}_stereoisomers.sdf"),
        sanitize=True,
        removeHs=False,
    )
    assert sum(1 for mol in supplier if mol is not None) == 3


def test_enumeration_ceiling_is_hard_bounded() -> None:
    config = RunConfig(stereoisomer_filtering={"enumeration_ceiling_multiplier": 1000})
    assert _enumeration_ceiling(config, 16) == STEREO_ENUMERATION_HARD_CEILING
    small = RunConfig(stereoisomer_filtering={"enumeration_ceiling_multiplier": 2})
    assert _enumeration_ceiling(small, 16) == 32
    with pytest.raises(ValueError, match="positive"):
        RunConfig(stereoisomer_filtering={"enumeration_ceiling_multiplier": 0})


def _stereo(
    smiles: str,
    tautomer_id: str,
    stereo_index: int,
    *,
    selected: bool = True,
) -> StereoRecord:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    metadata = {} if selected else {"selected": False, "rejection_reason": "beyond cap"}
    return StereoRecord(
        id=f"{tautomer_id}_s{stereo_index:02d}_x",
        parent_id=tautomer_id,
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        formal_charge=Chem.GetFormalCharge(mol),
        source_software="test",
        source_python_function="test",
        stereo_index=stereo_index,
        metadata=metadata,
        rdkit_mol=mol,
    )


def test_cross_tautomer_duplicate_guard_dedupes_and_refills(tmp_path: Path, monkeypatch) -> None:
    # Tautomer t01 enumerates A+B; tautomer t02 enumerates a copy of A plus C;
    # t02's unused pool holds D.
    candidates = [
        _stereo("C[C@H](O)C(=O)O", "mol_p01_t01", 1),
        _stereo("F[C@H](Cl)Br", "mol_p01_t01", 2),
        _stereo("C[C@@H](O)C(=O)O", "mol_p01_t02", 1),  # different record, same... no:
        _stereo("C[C@H](O)C(=O)O", "mol_p01_t02", 2),
    ]
    # Make exactly the t02 copy identical to the t01 A record's structure.
    # (Enantiomer of A at index 1 keeps keys distinct; the index-2 record is
    # the exact duplicate.)
    pool = [_stereo("C[C@H](O)C(=O)N", "mol_p01_t02", 3, selected=False)]
    offered_to_auto3d: list[str] = []

    def fake_run_auto3d(input_path, output_dir, **kwargs):
        ids = [line.split()[1] for line in input_path.read_text(encoding="utf-8").splitlines()]
        offered_to_auto3d.extend(ids)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = output_dir / "auto3d_output.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for line in input_path.read_text(encoding="utf-8").splitlines():
            smiles, line_id = line.split(maxsplit=1)
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", line_id)
            mol.SetProp("E_kcal_mol", "0.0")
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "mock"]

    monkeypatch.setattr(stereo_auto3d_filter, "run_auto3d", fake_run_auto3d)
    config = RunConfig(
        input_path=tmp_path / "in.smi",
        output_dir=tmp_path / "run",
        enumeration={"max_stereoisomers_per_tautomer": 2},
        stereoisomer_filtering={"max_stereoisomers_per_tautomer": 2},
    )

    result = filter_stereoisomers_with_auto3d([*candidates, *pool], config)

    # t02's exact duplicate of t01's "C[C@H](O)C(=O)O" is eliminated pre-ranking
    # and refilled from t02's own unused pool ("C[C@H](O)C(=O)N").
    duplicate = candidates[3]
    promoted = pool[0]
    rejected_by_id = {record.id: record for record in result.rejected_records}
    assert duplicate.id in rejected_by_id
    decision_reason = rejected_by_id[duplicate.id].metadata["stereo_energy_filtering"]["reason"]
    assert "rejected_exact_duplicate_pre_ranking" in decision_reason
    assert duplicate.id not in offered_to_auto3d
    promoted_record = next(record for record in result.all_records if record.id == promoted.id)
    assert promoted_record.metadata["selection_reason"] == "selected_by_stereo_dedupe_refill"
    assert promoted_record.metadata["stereo_refill"]["promoted_after_cross_tautomer_dedupe"] is True
    assert promoted.id in offered_to_auto3d
    # Only one copy of each exact structure is selected downstream.
    selected_smiles = [record.isomeric_smiles for record in result.selected_records]
    assert len(selected_smiles) == len(set(selected_smiles))

    audit_path = tmp_path / "run" / "stereoisomer_filtering" / "stereo_dedupe.csv"
    with audit_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    merge_rows = [row for row in rows if row["action"] == "merge"]
    refill_rows = [row for row in rows if row["action"] == "refill"]
    assert len(merge_rows) == 1
    assert merge_rows[0]["eliminated_stereo_ids"] == duplicate.id
    assert merge_rows[0]["retained_stereo_id"] == candidates[0].id
    assert len(refill_rows) == 1
    assert refill_rows[0]["promoted_stereo_id"] == promoted.id
    assert refill_rows[0]["parent_tautomer_id"] == "mol_p01_t02"
