from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field


class InputMolecule(BaseModel):
    index: int
    name: str
    input_format: str
    smiles: str | None = None
    source_path: Path
    input_hash: str
    properties: dict[str, str] = Field(default_factory=dict)


class ToolStatus(BaseModel):
    name: str
    kind: str
    required: bool = False
    available: bool
    detail: str = ""


class VariantRecord(BaseModel):
    variant_id: str
    parent_name: str
    smiles: str | None = None
    relative_energy_kcal_mol: float | None = None
    approximate_population: float | None = None
    status: str = "placeholder"


class WorkflowResult(BaseModel):
    outdir: Path
    molecule_count: int
    ranked_records: list[VariantRecord]

