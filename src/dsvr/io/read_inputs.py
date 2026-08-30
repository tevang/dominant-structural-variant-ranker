from __future__ import annotations

import csv
from pathlib import Path
from typing import Literal

from dsvr.io.sdf import read_sdf
from dsvr.io.smiles import read_smiles
from dsvr.models import MoleculeInput

InputFormat = Literal["auto", "smi", "smiles", "sdf"]

SMILES_SUFFIXES = {".smi", ".smiles", ".txt", ".csv"}
SDF_SUFFIXES = {".sdf", ".sd"}

INVALID_INPUT_COLUMNS = [
    "input_id",
    "source_format",
    "line_number",
    "name",
    "smiles",
    "raw_record",
    "error",
]


def read_molecules(
    path: Path,
    *,
    input_format: InputFormat = "auto",
    deduplicate: bool = True,
    invalid_output_path: Path | None = None,
) -> list[MoleculeInput]:
    molecules, _invalid_records = validate_input_file(
        path,
        input_format=input_format,
        deduplicate=deduplicate,
        invalid_output_path=invalid_output_path,
    )
    return molecules


def validate_input_file(
    path: Path,
    *,
    input_format: InputFormat = "auto",
    deduplicate: bool = True,
    invalid_output_path: Path | None = None,
) -> tuple[list[MoleculeInput], list[dict[str, str]]]:
    molecules, invalid_records = _read_by_format(
        path,
        input_format=input_format,
        deduplicate=deduplicate,
    )
    write_invalid_inputs_csv(
        invalid_output_path or path.parent / "invalid_inputs.csv",
        invalid_records,
    )
    return molecules, invalid_records


def _read_by_format(
    path: Path,
    *,
    input_format: InputFormat,
    deduplicate: bool,
) -> tuple[list[MoleculeInput], list[dict[str, str]]]:
    source_format = resolve_input_format(path, input_format)
    if source_format == "sdf":
        return read_sdf(path, deduplicate=deduplicate)
    return read_smiles(path, deduplicate=deduplicate)


def resolve_input_format(
    path: Path,
    input_format: InputFormat = "auto",
) -> Literal["smiles", "sdf"]:
    normalized = input_format.lower()
    if normalized in {"smi", "smiles"}:
        return "smiles"
    if normalized == "sdf":
        return "sdf"
    if normalized != "auto":
        raise ValueError(f"Unsupported input format: {input_format}")

    suffix = path.suffix.lower()
    if suffix in SMILES_SUFFIXES:
        return "smiles"
    if suffix in SDF_SUFFIXES:
        return "sdf"
    raise ValueError(f"Unsupported input extension: {path}")


def write_invalid_inputs_csv(path: Path, invalid_records: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INVALID_INPUT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(invalid_records)
