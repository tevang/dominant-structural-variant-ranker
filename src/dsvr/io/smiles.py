from __future__ import annotations

from pathlib import Path

from rdkit import Chem

from dsvr.models import MoleculeInput
from dsvr.utils.hashing import sha256_text

HEADER_SMILES_NAMES = {"smiles", "smile"}
HEADER_NAME_NAMES = {"molname", "name", "molecule", "molecule_name", "id"}


def read_smiles(
    path: Path,
    *,
    deduplicate: bool = True,
) -> tuple[list[MoleculeInput], list[dict[str, str]]]:
    molecules: list[MoleculeInput] = []
    invalid_records: list[dict[str, str]] = []
    seen_isomeric_smiles: set[str] = set()
    header: list[str] | None = None
    smiles_column = 0
    name_column: int | None = None
    data_index = 0

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if header is None and _is_header(fields):
                header = [field.lower() for field in fields]
                smiles_column = _find_column(header, HEADER_SMILES_NAMES, default=0)
                name_column = _find_column(header, HEADER_NAME_NAMES, default=None)
                continue

            data_index += 1
            input_id = f"mol_{data_index:06d}"
            try:
                original_smiles = fields[smiles_column]
            except IndexError:
                invalid_records.append(
                    _invalid_record(input_id, line_number, line, "missing SMILES column")
                )
                continue
            original_name = _extract_name(fields, smiles_column, name_column)
            molname = original_name or input_id
            molecule = Chem.MolFromSmiles(original_smiles)
            if molecule is None:
                invalid_records.append(
                    _invalid_record(input_id, line_number, line, "RDKit failed to parse SMILES")
                )
                continue

            canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
            isomeric_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
            warnings: list[str] = []
            if deduplicate and isomeric_smiles in seen_isomeric_smiles:
                warnings.append("duplicate canonical isomeric SMILES; skipped")
                continue
            seen_isomeric_smiles.add(isomeric_smiles)
            molecules.append(
                MoleculeInput(
                    input_id=input_id,
                    molname=molname,
                    source_format="smiles",
                    original_smiles=original_smiles,
                    canonical_smiles=canonical_smiles,
                    isomeric_smiles=isomeric_smiles,
                    rdkit_mol=molecule,
                    input_properties={
                        "line_number": str(line_number),
                        "record_index": str(len(molecules)),
                        "original_name": original_name,
                        "raw_record": line,
                        "input_hash": sha256_text(line),
                    },
                    warnings=warnings,
                )
            )
    return molecules, invalid_records


def _is_header(fields: list[str]) -> bool:
    lowered = [field.lower() for field in fields]
    return bool(lowered) and lowered[0] in HEADER_SMILES_NAMES


def _find_column(header: list[str], names: set[str], default: int | None) -> int | None:
    for index, name in enumerate(header):
        if name in names:
            return index
    return default


def _extract_name(fields: list[str], smiles_column: int, name_column: int | None) -> str:
    if name_column is not None:
        return fields[name_column].strip() if name_column < len(fields) else ""
    if smiles_column == 0 and len(fields) > 1:
        return " ".join(fields[1:]).strip()
    return ""


def _invalid_record(input_id: str, line_number: int, raw_record: str, error: str) -> dict[str, str]:
    return {
        "input_id": input_id,
        "source_format": "smiles",
        "line_number": str(line_number),
        "raw_record": raw_record,
        "error": error,
    }
