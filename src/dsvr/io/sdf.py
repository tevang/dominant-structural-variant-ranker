from __future__ import annotations

from pathlib import Path

from rdkit import Chem

from dsvr.models import MoleculeInput
from dsvr.utils.hashing import sha256_text


def read_sdf(
    path: Path,
    *,
    deduplicate: bool = True,
) -> tuple[list[MoleculeInput], list[dict[str, str]]]:
    records = _split_sdf_records(path)
    supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
    molecules: list[MoleculeInput] = []
    invalid_records: list[dict[str, str]] = []
    seen_isomeric_smiles: set[str] = set()

    for record_number, record in enumerate(records, start=1):
        input_id = f"mol_{record_number:06d}"
        molecule = supplier[record_number - 1] if record_number - 1 < len(supplier) else None
        if molecule is None:
            invalid_records.append(
                {
                    "input_id": input_id,
                    "source_format": "sdf",
                    "line_number": "",
                    "raw_record": record,
                    "error": "RDKit failed to parse or sanitize SDF record",
                }
            )
            continue

        properties = {
            name: molecule.GetProp(name)
            for name in molecule.GetPropNames(includePrivate=True)
        }
        raw_name = molecule.GetProp("_Name").strip() if molecule.HasProp("_Name") else ""
        molname = raw_name or input_id
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        warnings: list[str] = []
        if deduplicate and isomeric_smiles in seen_isomeric_smiles:
            warnings.append("duplicate canonical isomeric SMILES; skipped")
            continue
        seen_isomeric_smiles.add(isomeric_smiles)
        properties.update(
            {
                "record_index": str(len(molecules)),
                "record_number": str(record_number),
                "original_name": raw_name,
                "raw_record": record,
                "input_hash": sha256_text(record),
            }
        )
        molecules.append(
            MoleculeInput(
                input_id=input_id,
                molname=molname,
                source_format="sdf",
                original_smiles=None,
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                rdkit_mol=molecule,
                input_properties=properties,
                warnings=warnings,
            )
        )
    return molecules, invalid_records


def _split_sdf_records(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    return [record.strip() for record in text.split("$$$$") if record.strip()]
