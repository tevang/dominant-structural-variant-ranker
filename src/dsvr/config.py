from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel, Field, field_validator


class WorkflowConfig(BaseModel):
    ph: float = Field(default=7.0, description="Candidate-generation pH.")
    solvent: str = "water"
    solvent_model: str = "alpb"
    max_tautomers: int = 64
    max_stereoisomers: int = 64
    max_conformers: int = 50
    conformer_backend: str = "rdkit"
    auto3d_internal_enumeration: bool = False
    execute_external: bool = False
    optional_refinement: str = "none"

    @field_validator("max_tautomers", "max_stereoisomers", "max_conformers")
    @classmethod
    def positive_int(cls, value: int) -> int:
        if value < 1:
            raise ValueError("limits must be positive")
        return value


class DsvrConfig(BaseModel):
    workflow: WorkflowConfig = Field(default_factory=WorkflowConfig)


def load_config(path: Path) -> DsvrConfig:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return DsvrConfig.model_validate(data)

