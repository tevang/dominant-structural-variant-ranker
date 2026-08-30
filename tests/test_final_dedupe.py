"""Tests for the final safety-net exact-duplicate removal on 3D variants."""

from __future__ import annotations

import csv
from pathlib import Path

from rdkit import Chem

from dsvr.chemistry import final3d
from dsvr.chemistry.final3d import generate_final_3d_variants
from dsvr.config import RunConfig
from dsvr.models import StereoRecord


def _stereo(record_id: str, tautomer_id: str, smiles: str) -> StereoRecord:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return StereoRecord(
        id=record_id,
        parent_id=tautomer_id,
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        formal_charge=Chem.GetFormalCharge(mol),
        source_software="test",
        source_python_function="test",
        stereo_index=1,
        rdkit_mol=mol,
    )


def test_final_variants_dedupe_duplicate_geometry_records(tmp_path: Path, monkeypatch) -> None:
    energies = {
        "mol_p01_t01_s01_aaa": -5.0,
        "mol_p02_t01_s01_bbb": -10.0,
        "mol_p01_t01_s02_ccc": 0.0,
    }

    def fake_run_auto3d(input_path, output_dir, **kwargs):
        output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = output_dir / "mock_final.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        supplier = Chem.SDMolSupplier(str(input_path), sanitize=True, removeHs=False)
        for mol in supplier:
            if mol is None:
                continue
            stereo_id = mol.GetProp("DSVR_STEREO_ID")
            mol.SetProp("E_kcal_mol", str(energies[stereo_id]))
            # Give every duplicate-branch structure the exact same 3D geometry,
            # like converged Auto3D outputs from different protomer branches.
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "mock"]

    monkeypatch.setattr(final3d, "run_auto3d", fake_run_auto3d)
    stereo_records = [
        _stereo("mol_p01_t01_s01_aaa", "mol_p01_t01", "CCO"),
        _stereo("mol_p02_t01_s01_bbb", "mol_p02_t01", "CCO"),
        _stereo("mol_p01_t01_s02_ccc", "mol_p01_t01", "CCN"),
    ]
    config = RunConfig(input_path=tmp_path / "in.smi", output_dir=tmp_path / "run")

    result = generate_final_3d_variants(stereo_records, config)

    ids = {record.parent_id: record for record in result.records}
    # One output per unique structure: "CCO" appears once (lowest energy wins),
    # "CCN" stays.
    assert set(ids) == {"mol_p02_t01_s01_bbb", "mol_p01_t01_s02_ccc"}
    kept = ids["mol_p02_t01_s01_bbb"]
    assert kept.energy_kcal_mol == -10.0
    merged_from = kept.metadata["merged_from"]
    assert len(merged_from) == 1
    assert merged_from[0]["final_variant_id"].startswith("mol_p01_t01_s01_aaa")
    assert merged_from[0]["stereo_id"] == "mol_p01_t01_s01_aaa"
    assert merged_from[0]["input_molecule_id"] == "mol"

    audit_path = tmp_path / "run" / "final_dedupe_audit.csv"
    assert audit_path.exists()
    with audit_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["action"] == "merge"
    assert rows[0]["eliminated_final_variant_ids"].startswith("mol_p01_t01_s01_aaa")
    assert rows[0]["retained_final_variant_id"].startswith("mol_p02_t01_s01_bbb")

    # The final SDF mirrors the deduplicated records.
    final_sdf = tmp_path / "run" / "final_variants.sdf"
    supplier = Chem.SDMolSupplier(str(final_sdf), sanitize=True, removeHs=False)
    written = [mol for mol in supplier if mol is not None]
    assert len(written) == 2


def test_final_dedupe_can_be_disabled(tmp_path: Path, monkeypatch) -> None:
    def fake_run_auto3d(input_path, output_dir, **kwargs):
        output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = output_dir / "mock_final.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        supplier = Chem.SDMolSupplier(str(input_path), sanitize=True, removeHs=False)
        for mol in supplier:
            if mol is None:
                continue
            mol.SetProp("E_kcal_mol", "0.0")
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "mock"]

    monkeypatch.setattr(final3d, "run_auto3d", fake_run_auto3d)
    stereo_records = [
        _stereo("mol_p01_t01_s01_aaa", "mol_p01_t01", "CCO"),
        _stereo("mol_p02_t01_s01_bbb", "mol_p02_t01", "CCO"),
    ]
    config = RunConfig(
        input_path=tmp_path / "in.smi",
        output_dir=tmp_path / "run",
        final_3d={"dedupe_final_variants": False},
    )

    result = generate_final_3d_variants(stereo_records, config)

    assert len(result.records) == 2
    audit_path = tmp_path / "run" / "final_dedupe_audit.csv"
    with audit_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows == []
