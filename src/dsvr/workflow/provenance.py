from __future__ import annotations

import csv
import json
import platform
import sys
from pathlib import Path
from typing import Any, TypeVar

from pydantic import BaseModel

from dsvr.config import RunConfig
from dsvr.models import (
    AnyLineageRecord,
    CrestConformerRecord,
    MoleculeInput,
    MoleculeRecord,
    ProtomerRecord,
    RankedVariantRecord,
    SeedConformerRecord,
    StereoRecord,
    TautomerRecord,
    ThermoRecord,
)

LineageT = TypeVar("LineageT", bound=AnyLineageRecord)


def build_provenance(
    input_path: Path,
    config: RunConfig,
    molecules: list[MoleculeInput],
) -> dict[str, object]:
    return {
        "input_path": str(input_path),
        "config": config.model_dump(mode="json"),
        "molecule_count": len(molecules),
        "input_hashes": [mol.input_hash for mol in molecules],
        "python": sys.version,
        "platform": platform.platform(),
        "scientific_limitation": (
            "pH controls candidate generation by default; cross-protomer "
            "populations are approximate "
            "without explicit micro-pKa/proton chemical-potential corrections."
        ),
    }


def write_provenance_jsonl(records: list[AnyLineageRecord], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_jsonl(
        output_dir / "enumeration_provenance.jsonl",
        [
            record
            for record in records
            if isinstance(record, (MoleculeRecord, ProtomerRecord, TautomerRecord, StereoRecord))
        ],
    )
    write_jsonl(
        output_dir / "conformer_provenance.jsonl",
        [
            record
            for record in records
            if isinstance(record, (SeedConformerRecord, CrestConformerRecord, ThermoRecord))
        ],
    )
    write_jsonl(
        output_dir / "ranking_provenance.jsonl",
        [record for record in records if isinstance(record, RankedVariantRecord)],
    )


def write_summary_tables(records: list[AnyLineageRecord], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "inputs.csv", _records_of_type(records, MoleculeRecord))
    _write_csv(output_dir / "protomers.csv", _records_of_type(records, ProtomerRecord))
    _write_csv(output_dir / "tautomers.csv", _records_of_type(records, TautomerRecord))
    _write_csv(output_dir / "stereoisomers.csv", _records_of_type(records, StereoRecord))
    _write_csv(output_dir / "seeds.csv", _records_of_type(records, SeedConformerRecord))
    _write_csv(output_dir / "crest_conformers.csv", _records_of_type(records, CrestConformerRecord))
    _write_csv(output_dir / "thermo.csv", _records_of_type(records, ThermoRecord))
    _write_csv(output_dir / "ranked_variants.csv", _records_of_type(records, RankedVariantRecord))


def write_all_provenance_outputs(records: list[AnyLineageRecord], output_dir: Path) -> None:
    write_provenance_jsonl(records, output_dir)
    write_summary_tables(records, output_dir)


def write_jsonl(path: Path, records: list[AnyLineageRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(_jsonable(record), sort_keys=True) + "\n")


def _write_csv(path: Path, records: list[AnyLineageRecord]) -> None:
    columns = _columns_for_records(records)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            row = _flatten_record(record)
            writer.writerow(row)


def _records_of_type(
    records: list[AnyLineageRecord],
    record_type: type[LineageT],
) -> list[LineageT]:
    return [record for record in records if isinstance(record, record_type)]


def _columns_for_records(records: list[AnyLineageRecord]) -> list[str]:
    base_columns = [
        "id",
        "parent_id",
        "input_molecule_id",
        "molname",
        "canonical_smiles",
        "isomeric_smiles",
        "molecular_formula",
        "formal_charge",
        "explicit_proton_count",
        "stage_name",
        "source_software",
        "source_command",
        "source_python_function",
        "output_paths",
        "warnings",
        "metadata",
    ]
    extra_columns: list[str] = []
    for record in records:
        for key in record.model_dump(mode="json"):
            if key not in base_columns and key not in extra_columns:
                extra_columns.append(key)
    return [*base_columns, *extra_columns]


def _flatten_record(record: AnyLineageRecord) -> dict[str, Any]:
    data = record.model_dump(mode="json")
    for key in ("output_paths", "warnings"):
        data[key] = json.dumps(data.get(key, []), sort_keys=True)
    data["metadata"] = json.dumps(data.get("metadata", {}), sort_keys=True)
    return data


def _jsonable(value: Any) -> Any:
    if isinstance(value, BaseModel):
        return value.model_dump(mode="json")
    if isinstance(value, Path):
        return str(value)
    return value
