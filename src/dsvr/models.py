from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

from dsvr.utils.hashing import sha256_text

StageName = Literal[
    "input",
    "protomer",
    "tautomer",
    "stereo",
    "seed_conformer",
    "crest_conformer",
    "thermo",
    "ranked_variant",
]


class MoleculeInput(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    input_id: str
    molname: str
    source_format: str
    original_smiles: str | None
    canonical_smiles: str
    isomeric_smiles: str
    rdkit_mol: Any
    input_properties: dict[str, str] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)

    @property
    def index(self) -> int:
        raw_index = self.input_properties.get("record_index")
        return int(raw_index) if raw_index is not None else 0

    @property
    def name(self) -> str:
        return self.molname

    @property
    def input_format(self) -> str:
        return self.source_format

    @property
    def smiles(self) -> str:
        return self.original_smiles or self.isomeric_smiles

    @property
    def input_hash(self) -> str:
        return self.input_id

    @property
    def properties(self) -> dict[str, str]:
        return self.input_properties


class LineageRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str
    parent_id: str | None
    input_molecule_id: str
    molname: str
    canonical_smiles: str | None = None
    isomeric_smiles: str | None = None
    molecular_formula: str | None = None
    formal_charge: int | None = None
    explicit_proton_count: int | None = None
    stage_name: StageName
    source_software: str
    source_command: str | None = None
    source_python_function: str | None = None
    output_paths: list[Path] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)

    @field_validator("output_paths", mode="before")
    @classmethod
    def normalize_output_paths(cls, value: object) -> object:
        if value is None:
            return []
        if isinstance(value, (str, Path)):
            return [value]
        return value


class MoleculeRecord(LineageRecord):
    stage_name: Literal["input"] = "input"


class ProtomerRecord(LineageRecord):
    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    stage_name: Literal["protomer"] = "protomer"
    protomer_index: int = 0
    rdkit_mol: Any = Field(default=None, exclude=True)


class TautomerRecord(LineageRecord):
    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    stage_name: Literal["tautomer"] = "tautomer"
    tautomer_index: int = 0
    rdkit_mol: Any = Field(default=None, exclude=True)


class StereoRecord(LineageRecord):
    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    stage_name: Literal["stereo"] = "stereo"
    stereo_index: int = 0
    rdkit_mol: Any = Field(default=None, exclude=True)


class SeedConformerRecord(LineageRecord):
    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    stage_name: Literal["seed_conformer"] = "seed_conformer"
    conformer_index: int = 0
    energy_kcal_mol: float | None = None
    rdkit_mol: Any = Field(default=None, exclude=True)
    rdkit_conformer_id: int | None = None
    forcefield: str | None = None
    forcefield_status: str = "not_run"
    minimization_converged: bool | None = None
    embedding_status: str = "not_run"


class CrestConformerRecord(LineageRecord):
    stage_name: Literal["crest_conformer"] = "crest_conformer"
    crest_index: int = 0
    energy_kcal_mol: float | None = None
    relative_energy_kcal_mol: float | None = None


class ThermoRecord(LineageRecord):
    stage_name: Literal["thermo"] = "thermo"
    temperature_kelvin: float = 298.15
    free_energy_kcal_mol: float | None = None
    entropy_cal_mol_k: float | None = None


class RankedVariantRecord(LineageRecord):
    stage_name: Literal["ranked_variant"] = "ranked_variant"
    rank: int
    score_kcal_mol: float | None = None
    relative_free_energy_kcal_mol: float | None = None
    boltzmann_population: float | None = None
    population_scope: str = "same_formula"
    approximate_population: bool = True


AnyLineageRecord = (
    MoleculeRecord
    | ProtomerRecord
    | TautomerRecord
    | StereoRecord
    | SeedConformerRecord
    | CrestConformerRecord
    | ThermoRecord
    | RankedVariantRecord
)


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
    version: str | None = None
    minimum_version: str | None = None
    meets_minimum_version: bool | None = None


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


def make_input_id(
    molname: str,
    canonical_smiles: str | None,
    metadata: dict[str, Any] | None = None,
) -> str:
    digest = short_hash(_hash_payload(canonical_smiles, None, metadata))
    return f"{sanitize_id_part(molname)}_{digest}"


def make_protomer_id(
    input_id: str,
    protomer_index: int,
    canonical_smiles: str | None,
    isomeric_smiles: str | None,
    metadata: dict[str, Any] | None = None,
) -> str:
    return _child_id(
        input_id,
        f"p{protomer_index:02d}",
        canonical_smiles,
        isomeric_smiles,
        metadata,
    )


def make_tautomer_id(
    protomer_id: str,
    tautomer_index: int,
    canonical_smiles: str | None,
    isomeric_smiles: str | None,
    metadata: dict[str, Any] | None = None,
) -> str:
    return _child_id(
        protomer_id,
        f"t{tautomer_index:02d}",
        canonical_smiles,
        isomeric_smiles,
        metadata,
    )


def make_stereo_id(
    tautomer_id: str,
    stereo_index: int,
    canonical_smiles: str | None,
    isomeric_smiles: str | None,
    metadata: dict[str, Any] | None = None,
) -> str:
    return _child_id(
        tautomer_id,
        f"s{stereo_index:02d}",
        canonical_smiles,
        isomeric_smiles,
        metadata,
    )


def make_seed_id(
    stereo_id: str,
    conformer_index: int,
    canonical_smiles: str | None,
    isomeric_smiles: str | None,
    metadata: dict[str, Any] | None = None,
) -> str:
    return _child_id(
        stereo_id,
        f"c{conformer_index:02d}",
        canonical_smiles,
        isomeric_smiles,
        metadata,
    )


def make_crest_conformer_id(
    stereo_id: str,
    crest_index: int,
    canonical_smiles: str | None,
    isomeric_smiles: str | None,
    metadata: dict[str, Any] | None = None,
) -> str:
    return _child_id(
        stereo_id,
        f"crest{crest_index:02d}",
        canonical_smiles,
        isomeric_smiles,
        metadata,
    )


def make_thermo_id(
    conformer_id: str,
    canonical_smiles: str | None,
    isomeric_smiles: str | None,
    metadata: dict[str, Any] | None = None,
) -> str:
    return _child_id(conformer_id, "thermo", canonical_smiles, isomeric_smiles, metadata)


def make_ranked_variant_id(
    parent_id: str,
    rank: int,
    canonical_smiles: str | None,
    isomeric_smiles: str | None,
    metadata: dict[str, Any] | None = None,
) -> str:
    return _child_id(parent_id, f"rank{rank:04d}", canonical_smiles, isomeric_smiles, metadata)


def sanitize_id_part(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_").lower()
    return cleaned or "mol"


def short_hash(value: str, length: int = 10) -> str:
    return sha256_text(value)[:length]


def _child_id(
    parent_id: str,
    tag: str,
    canonical_smiles: str | None,
    isomeric_smiles: str | None,
    metadata: dict[str, Any] | None,
) -> str:
    digest = short_hash(_hash_payload(canonical_smiles, isomeric_smiles, metadata))
    return f"{parent_id}_{tag}_{digest}"


def _hash_payload(
    canonical_smiles: str | None,
    isomeric_smiles: str | None,
    metadata: dict[str, Any] | None,
) -> str:
    payload = {
        "canonical_smiles": canonical_smiles or "",
        "isomeric_smiles": isomeric_smiles or "",
        "metadata": metadata or {},
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
