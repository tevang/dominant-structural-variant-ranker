from __future__ import annotations

from collections import Counter
from pathlib import Path

from rdkit import Chem

from dsvr.models import MoleculeInput
from dsvr.utils.hashing import sha256_text

HEADER_SMILES_NAMES = {
    "smiles",
    "smile",
    "canonical_smiles",
    "isomeric_smiles",
    "smiles_string",
    "structure_smiles",
}
HEADER_NAME_NAMES = {
    "molname",
    "name",
    "molecule",
    "molecule_name",
    "id",
    "chembl_id",
    "molecule chembl id",
    "molecule_chembl_id",
    "cmpd_id",
    "compound_id",
}

NO_SMILES_COLUMN_ERROR = "no SMILES column found"

_DELIMITER_CANDIDATES: tuple[str | None, ...] = (",", ";", "\t", None)
_DELIMITER_SAMPLE_LINES = 50
_PROBE_SAMPLE_ROWS = 25


def read_smiles(
    path: Path,
    *,
    deduplicate: bool = True,
) -> tuple[list[MoleculeInput], list[dict[str, str]]]:
    molecules: list[MoleculeInput] = []
    invalid_records: list[dict[str, str]] = []
    seen_isomeric_smiles: set[str] = set()
    rows = _data_rows(path)
    if not rows:
        return molecules, invalid_records

    delimiter = _sniff_delimiter([line for _, line in rows])
    header, smiles_column, name_column = _resolve_columns(rows, delimiter)
    if header is not None:
        rows = rows[1:]

    for data_index, (line_number, line) in enumerate(rows, start=1):
        input_id = f"mol_{data_index:06d}"
        fields = _split_fields(line, delimiter)
        if smiles_column is None:
            invalid_records.append(
                _invalid_record(
                    input_id,
                    line_number,
                    line,
                    name=_extract_name(fields, None, name_column),
                    smiles="",
                    error=NO_SMILES_COLUMN_ERROR,
                )
            )
            continue
        original_name = _extract_name(fields, smiles_column, name_column)
        if smiles_column >= len(fields):
            invalid_records.append(
                _invalid_record(
                    input_id,
                    line_number,
                    line,
                    name=original_name,
                    smiles="",
                    error="missing SMILES field",
                )
            )
            continue
        original_smiles = fields[smiles_column]
        molname = original_name or input_id
        molecule = Chem.MolFromSmiles(original_smiles)
        if molecule is None:
            invalid_records.append(
                _invalid_record(
                    input_id,
                    line_number,
                    line,
                    name=original_name,
                    smiles=original_smiles,
                    error="RDKit failed to parse SMILES",
                )
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


def _data_rows(path: Path) -> list[tuple[int, str]]:
    rows: list[tuple[int, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append((line_number, line))
    return rows


def _split_fields(line: str, delimiter: str | None) -> list[str]:
    if delimiter is None:
        return line.split()
    return [field.strip() for field in line.split(delimiter)]


def _sniff_delimiter(lines: list[str]) -> str | None:
    sample = lines[:_DELIMITER_SAMPLE_LINES]
    best: str | None = None
    best_score = -1.0
    for candidate in _DELIMITER_CANDIDATES:
        counts = [len(_split_fields(line, candidate)) for line in sample]
        mode_count, mode = Counter(counts).most_common(1)[0]
        if candidate is not None and mode < 2:
            continue
        score = mode_count / len(counts)
        if score > best_score:
            best, best_score = candidate, score
    return best


def _resolve_columns(
    rows: list[tuple[int, str]],
    delimiter: str | None,
) -> tuple[list[str] | None, int | None, int | None]:
    header: list[str] | None = None
    smiles_column: int | None = None
    name_column: int | None = None
    raw_first_fields = _split_fields(rows[0][1], delimiter)
    first_fields = [field.strip().lower() for field in raw_first_fields]
    if _looks_like_header(first_fields, raw_first_fields):
        header = first_fields
        smiles_column = _find_column(header, HEADER_SMILES_NAMES, default=None)
        name_column = _find_column(header, HEADER_NAME_NAMES, default=None)

    if smiles_column is None:
        data_rows = rows[1:] if header is not None else rows
        excluded = (
            {index for index, cell in enumerate(header) if cell in HEADER_NAME_NAMES}
            if header is not None
            else set()
        )
        smiles_column = _probe_smiles_column(data_rows, delimiter, excluded)
    return header, smiles_column, name_column


def _looks_like_header(normalized: list[str], original: list[str]) -> bool:
    if any(cell in HEADER_SMILES_NAMES for cell in normalized):
        return True
    if not any(cell in HEADER_NAME_NAMES for cell in normalized):
        return False
    # A line with a recognized name cell is still data when it also carries a
    # parseable SMILES (e.g. a molecule actually named "id"). Probe original
    # case: header matching is case-insensitive, SMILES parsing is not.
    return not any(_is_parseable_smiles(cell) for cell in original)


def _probe_smiles_column(
    rows: list[tuple[int, str]],
    delimiter: str | None,
    excluded: set[int],
) -> int | None:
    sample = [_split_fields(line, delimiter) for _, line in rows[:_PROBE_SAMPLE_ROWS]]
    column_count = max((len(fields) for fields in sample), default=0)
    best_column: int | None = None
    best_hits = 0
    for column in range(column_count):
        if column in excluded:
            continue
        hits = sum(
            1
            for fields in sample
            if column < len(fields) and _is_parseable_smiles(fields[column])
        )
        if hits > best_hits:
            best_column, best_hits = column, hits
    return best_column


def _is_parseable_smiles(text: str) -> bool:
    return bool(text) and Chem.MolFromSmiles(text) is not None


def _find_column(header: list[str], names: set[str], default: int | None) -> int | None:
    for index, name in enumerate(header):
        if name in names:
            return index
    return default


def _extract_name(fields: list[str], smiles_column: int | None, name_column: int | None) -> str:
    if (
        name_column is not None
        and name_column != smiles_column
        and name_column < len(fields)
    ):
        return fields[name_column].strip()
    if smiles_column is not None and len(fields) > 1:
        return " ".join(
            field for index, field in enumerate(fields) if index != smiles_column
        ).strip()
    return ""


def _invalid_record(
    input_id: str,
    line_number: int,
    raw_record: str,
    *,
    name: str,
    smiles: str,
    error: str,
) -> dict[str, str]:
    return {
        "input_id": input_id,
        "source_format": "smiles",
        "line_number": str(line_number),
        "name": name,
        "smiles": smiles,
        "raw_record": raw_record,
        "error": error,
    }
