"""Unit tests for engine-level cross-branch tautomer dedupation and refill."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.models import TautomerRecord
from dsvr.reporting.progress import ProgressRecorder
from dsvr.workflow.engine import dedupe_and_refill_tautomers


def _tautomer(
    record_id: str,
    parent_id: str,
    smiles: str,
    *,
    input_id: str = "mol",
    relative_energy: float | None = None,
    selected: bool = True,
) -> TautomerRecord:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    metadata: dict = {}
    if relative_energy is not None:
        metadata["auto3d_tautomer_filtering"] = {
            "selected": selected,
            "relative_energy_kcal_mol": relative_energy,
            "rank": None,
        }
    if not selected:
        metadata["selected"] = False
        metadata["rejection_reason"] = "beyond cap"
    return TautomerRecord(
        id=record_id,
        parent_id=parent_id,
        input_molecule_id=input_id,
        molname="mol",
        canonical_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        formal_charge=Chem.GetFormalCharge(mol),
        source_software="test",
        tautomer_index=1,
        metadata=metadata,
        rdkit_mol=mol,
    )


def _read_audit(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_cross_branch_duplicates_merge_and_branch_refills_to_tauto_k(tmp_path: Path) -> None:
    config = RunConfig(output_dir=tmp_path / "run", tautomer_filtering={"tauto_k": 2})
    selected = [
        _tautomer("mol_p01_t01", "mol_p01", "CC(=O)OC", relative_energy=1.0),
        _tautomer("mol_p01_t02", "mol_p01", "CCO", relative_energy=0.5),
        _tautomer("mol_p02_t01", "mol_p02", "CC(=O)OC", relative_energy=0.2),
        _tautomer("mol_p02_t02", "mol_p02", "CCC", relative_energy=0.9),
    ]
    pool = [
        # Same structure as the selected "CCC": must be skipped during refill.
        _tautomer("mol_p01_t90", "mol_p01", "CCC", relative_energy=1.5, selected=False),
        _tautomer("mol_p01_t91", "mol_p01", "CCN", relative_energy=2.0, selected=False),
    ]

    result = dedupe_and_refill_tautomers(selected, pool, config)

    ids = {record.id for record in result}
    # Representative keeps the best (lowest) relative energy across branches.
    assert ids == {"mol_p02_t01", "mol_p01_t02", "mol_p02_t02", "mol_p01_t91"}
    representative = next(record for record in result if record.id == "mol_p02_t01")
    merged_from = representative.metadata["merged_from"]
    assert [entry["tautomer_id"] for entry in merged_from] == ["mol_p01_t01"]
    assert merged_from[0]["protomer_id"] == "mol_p01"
    # Refilled branch p01 back to tauto_k with the unique pool candidate.
    promoted = next(record for record in result if record.id == "mol_p01_t91")
    filtering = promoted.metadata["auto3d_tautomer_filtering"]
    assert filtering["selected"] is True
    assert filtering["reason"] == "selected_by_tautomer_dedupe_refill"
    assert promoted.metadata["tautomer_refill"]["promoted_after_cross_branch_dedupe"] is True

    rows = _read_audit(tmp_path / "run" / "enumeration" / "tautomers" / "tautomer_dedupe.csv")
    merge_rows = [row for row in rows if row["action"] == "merge"]
    refill_rows = [row for row in rows if row["action"] == "refill"]
    assert len(merge_rows) == 1
    assert merge_rows[0]["retained_tautomer_id"] == "mol_p02_t01"
    assert merge_rows[0]["eliminated_tautomer_ids"] == "mol_p01_t01"
    assert len(refill_rows) == 1
    assert refill_rows[0]["promoted_tautomer_id"] == "mol_p01_t91"
    assert refill_rows[0]["parent_protomer_id"] == "mol_p01"


def test_exhausted_refill_pool_records_shortfall(tmp_path: Path) -> None:
    config = RunConfig(output_dir=tmp_path / "run", tautomer_filtering={"tauto_k": 2})
    selected = [
        _tautomer("mol_p01_t01", "mol_p01", "c1ccncc1", relative_energy=0.1),
        _tautomer("mol_p01_t02", "mol_p01", "CCO", relative_energy=0.5),
        _tautomer("mol_p02_t01", "mol_p02", "c1ccncc1", relative_energy=0.9),
        _tautomer("mol_p02_t02", "mol_p02", "CCCl", relative_energy=1.1),
    ]

    progress = ProgressRecorder(tmp_path / "run", terminal=False)
    result = dedupe_and_refill_tautomers(selected, [], config, progress)

    ids = {record.id for record in result}
    assert ids == {"mol_p01_t01", "mol_p01_t02", "mol_p02_t02"}

    rows = _read_audit(tmp_path / "run" / "enumeration" / "tautomers" / "tautomer_dedupe.csv")
    shortfall_rows = [row for row in rows if row["action"] == "shortfall"]
    assert len(shortfall_rows) == 1
    assert shortfall_rows[0]["parent_protomer_id"] == "mol_p02"
    assert "1 of 2" in shortfall_rows[0]["detail"]

    warnings_text = (tmp_path / "run" / "warnings.jsonl").read_text(encoding="utf-8")
    assert "tautomer refill shortfall" in warnings_text
    assert "mol_p02" in warnings_text


def test_resume_loaded_records_dedupe_by_deterministic_id_tiebreak(tmp_path: Path) -> None:
    """Resume-loaded records carry no rank metadata: the smallest record ID
    wins deterministically, and branches cannot refill without a pool."""

    config = RunConfig(output_dir=tmp_path / "run", tautomer_filtering={"tauto_k": 2})
    selected = [
        _tautomer("mol_p01_t02_zzz", "mol_p01", "CC(=O)OC"),
        _tautomer("mol_p02_t01_aaa", "mol_p02", "CC(=O)OC"),
        _tautomer("mol_p01_t03_bbb", "mol_p01", "CCO"),
    ]

    result = dedupe_and_refill_tautomers(selected, [], config)

    ids = [record.id for record in result]
    assert "mol_p01_t02_zzz" in ids  # smallest ID wins the tiebreak
    assert "mol_p02_t01_aaa" not in ids
    assert ids.count("mol_p01_t02_zzz") == 1
    representative = next(record for record in result if record.id == "mol_p01_t02_zzz")
    assert [entry["tautomer_id"] for entry in representative.metadata["merged_from"]] == [
        "mol_p02_t01_aaa"
    ]


def test_distinct_tautomers_are_never_merged(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    config = RunConfig(output_dir=tmp_path / "run", tautomer_filtering={"tauto_k": 2})
    selected = [
        _tautomer("mol_p01_t01", "mol_p01", "O=C(NC1CCCCC1)c1cccc(-c2n[nH]c3ccc(-c4nc[nH]n4)cc23)c1"),
        _tautomer("mol_p02_t01", "mol_p02", "O=C(NC1CCCCC1)c1cccc(-c2nc3ccc(-c4nc[nH]n4)cc3[nH]2)c1"),
    ]

    result = dedupe_and_refill_tautomers(selected, [], config)

    assert {record.id for record in result} == {"mol_p01_t01", "mol_p02_t01"}
    rows = _read_audit(tmp_path / "run" / "enumeration" / "tautomers" / "tautomer_dedupe.csv")
    assert [row for row in rows if row["action"] == "merge"] == []
