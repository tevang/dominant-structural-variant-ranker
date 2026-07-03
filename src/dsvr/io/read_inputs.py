from __future__ import annotations

import csv
from pathlib import Path
from typing import Literal

from dsvr.io.sdf import read_sdf
from dsvr.io.smiles import read_smiles
from dsvr.models import MoleculeInput

InputFormat = Literal["auto", "smi", "smiles", "sdf"]

SMILES_SUFFIXES = {".smi", ".smiles", ".txt"}
SDF_SUFFIXES = {".sdf", ".sd"}


def read_molecules(
    path: Path,
    *,
    input_format: InputFormat = "auto",
    deduplicate: bool = True,
    invalid_output_path: Path | None = None,
) -> list[MoleculeInput]:
    source_format = resolve_input_format(path, input_format)
    if source_format == "sdf":
        molecules, invalid_records = read_sdf(path, deduplicate=deduplicate)
    else:
        molecules, invalid_records = read_smiles(path, deduplicate=deduplicate)

    if invalid_records:
        write_invalid_inputs_csv(
            invalid_output_path or path.parent / "invalid_inputs.csv",
            invalid_records,
        )
    return molecules


def validate_input_file(
    path: Path,
    *,
    input_format: InputFormat = "auto",
    deduplicate: bool = True,
    invalid_output_path: Path | None = None,
) -> tuple[list[MoleculeInput], list[dict[str, str]]]:
    source_format = resolve_input_format(path, input_format)
    if source_format == "sdf":
        molecules, invalid_records = read_sdf(path, deduplicate=deduplicate)
    else:
        molecules, invalid_records = read_smiles(path, deduplicate=deduplicate)
    if invalid_records:
        write_invalid_inputs_csv(
            invalid_output_path or path.parent / "invalid_inputs.csv",
            invalid_records,
        )
    return molecules, invalid_records


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
    columns = ["input_id", "source_format", "line_number", "raw_record", "error"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(invalid_records)
