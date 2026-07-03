from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

InputFormat = Literal["auto", "smi", "smiles", "sdf"]
SolventModel = Literal["alpb", "gbsa", "none"]
SeederMethod = Literal["etkdg", "auto3d", "both"]
RdkitForcefield = Literal["uff", "mmff", "none"]
Auto3dModel = Literal["AIMNet2", "ANI2x", "ANI2xt", "auto"]
PopulationScope = Literal["same_formula", "same_charge", "all_approximate"]
EnergyUnit = Literal["kcal/mol"]
LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR"]
QmBackend = Literal["psi4", "pyscf", "none"]

KNOWN_SOLVENTS = {
    "acetone",
    "acetonitrile",
    "benzene",
    "chloroform",
    "dmf",
    "dmso",
    "ethanol",
    "ether",
    "hexane",
    "methanol",
    "thf",
    "toluene",
    "water",
}


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ChemistryConfig(StrictModel):
    ph: float = 7.0
    ph_low: float | None = None
    ph_high: float | None = None
    solvent: str = "water"
    solvent_model: SolventModel = "alpb"
    temperature_kelvin: float = 298.15
    standardize: bool = True
    keep_salts: bool = False
    largest_fragment_only: bool = True

    @model_validator(mode="after")
    def validate_ph_window_and_solvent(self) -> ChemistryConfig:
        if self.ph_low is None and self.ph_high is None:
            object.__setattr__(self, "ph_low", self.ph)
            object.__setattr__(self, "ph_high", self.ph)
        elif self.ph_low is None or self.ph_high is None:
            raise ValueError("ph_low and ph_high must both be set or both be null")
        elif self.ph_low > self.ph_high:
            raise ValueError("ph_low must be <= ph_high")

        solvent_key = self.solvent.strip().lower()
        if not solvent_key:
            raise ValueError("solvent must not be empty")
        if solvent_key not in KNOWN_SOLVENTS:
            warnings.warn(
                f"Solvent '{self.solvent}' is not in DSVR's conservative known-solvent list; "
                "it will still be passed through to configured external tools.",
                stacklevel=2,
            )
        return self

    @field_validator("temperature_kelvin")
    @classmethod
    def positive_temperature(cls, value: float) -> float:
        if value <= 0:
            raise ValueError("temperature_kelvin must be positive")
        return value


class EnumerationConfig(StrictModel):
    max_protomers_per_molecule: int = 32
    max_tautomers_per_protomer: int = 64
    max_stereoisomers_per_tautomer: int = 64
    stereo_try_embedding: bool = True
    stereo_only_unassigned: bool = True
    stereo_unique: bool = True
    stereo_random_seed: int = 61453
    fail_on_enumeration_cap: bool = False

    @field_validator(
        "max_protomers_per_molecule",
        "max_tautomers_per_protomer",
        "max_stereoisomers_per_tautomer",
    )
    @classmethod
    def positive_cap(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("enumeration caps must be positive")
        return value


class SeedingConfig(StrictModel):
    method: SeederMethod = "etkdg"
    rdkit_num_conformers: int = 20
    rdkit_prune_rms_thresh: float = 0.5
    rdkit_forcefield: RdkitForcefield = "uff"
    auto3d_k: int = 5
    auto3d_model: Auto3dModel = "AIMNet2"
    auto3d_internal_tautomer_stereo_enum: bool = False

    @field_validator("rdkit_num_conformers", "auto3d_k")
    @classmethod
    def positive_count(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("conformer counts must be positive")
        return value

    @field_validator("rdkit_prune_rms_thresh")
    @classmethod
    def nonnegative_rms(cls, value: float) -> float:
        if value < 0:
            raise ValueError("rdkit_prune_rms_thresh must be non-negative")
        return value


class CrestConfig(StrictModel):
    enabled: bool = True
    executable: str = "crest"
    xtb_executable: str = "xtb"
    gfn: int = 2
    ewin_kcal_mol: float = 6.0
    nproc: int = 4
    walltime_minutes: int | None = None
    command_template: str | None = None
    extra_args: list[str] = Field(default_factory=list)

    @field_validator("gfn", "nproc")
    @classmethod
    def positive_int(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("gfn and nproc must be positive")
        return value

    @field_validator("ewin_kcal_mol")
    @classmethod
    def positive_energy_window(cls, value: float) -> float:
        if value <= 0:
            raise ValueError("ewin_kcal_mol must be positive")
        return value

    @field_validator("walltime_minutes")
    @classmethod
    def optional_positive_walltime(cls, value: int | None) -> int | None:
        if value is not None and value <= 0:
            raise ValueError("walltime_minutes must be positive when set")
        return value


class ThermoConfig(StrictModel):
    enabled: bool = True
    xtb_hessian: bool = True
    xtb_thermo: bool = True
    rrho_cutoff: float = 100.0
    population_scope: PopulationScope = "same_formula"
    energy_unit: EnergyUnit = "kcal/mol"

    @field_validator("rrho_cutoff")
    @classmethod
    def nonnegative_rrho_cutoff(cls, value: float) -> float:
        if value < 0:
            raise ValueError("rrho_cutoff must be non-negative")
        return value


class OptionalRefinementConfig(StrictModel):
    censo_enabled: bool = False
    censo_executable: str = "censo"
    censo_command_template: str | None = None
    censo_extra_args: list[str] = Field(default_factory=list)
    qm_backend: QmBackend = "none"
    qm_method: str = "b3lyp"
    qm_basis: str = "def2-svp"
    qm_solvent: str | None = None
    qm_optimize: bool = False
    qm_single_point: bool = True
    psi4_enabled: bool = False
    pyscf_enabled: bool = False
    max_candidates_for_refinement: int = 5

    @field_validator("max_candidates_for_refinement")
    @classmethod
    def positive_refinement_count(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("max_candidates_for_refinement must be positive")
        return value


class LoggingConfig(StrictModel):
    level: LogLevel = "INFO"
    save_subprocess_logs: bool = True
    tail_subprocess_logs: bool = True


class RunConfig(StrictModel):
    run_name: str = "dsvr-run"
    input_path: Path = Path("examples/test_molecules_minimal.smi")
    input_format: InputFormat = "auto"
    output_dir: Path = Path("runs/dsvr")
    max_workers: int = 1
    dry_run: bool = False
    overwrite: bool = False
    resume: bool = True
    chemistry: ChemistryConfig = Field(default_factory=ChemistryConfig)
    enumeration: EnumerationConfig = Field(default_factory=EnumerationConfig)
    seeding: SeedingConfig = Field(default_factory=SeedingConfig)
    crest: CrestConfig = Field(default_factory=CrestConfig)
    thermo: ThermoConfig = Field(default_factory=ThermoConfig)
    refinement: OptionalRefinementConfig = Field(default_factory=OptionalRefinementConfig)
    logging: LoggingConfig = Field(default_factory=LoggingConfig)

    @field_validator("max_workers")
    @classmethod
    def positive_max_workers(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("max_workers must be positive")
        return value


def load_config(path: Path) -> RunConfig:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return RunConfig.model_validate(data)


def merge_cli_overrides(config: RunConfig, **overrides: Any) -> RunConfig:
    data = config.model_dump(mode="python")
    _set_if_present(data, "input_path", overrides.get("input_path"))
    _set_if_present(data, "output_dir", overrides.get("output_dir"))
    _set_if_present(data["chemistry"], "ph", overrides.get("ph"))
    if overrides.get("ph") is not None:
        data["chemistry"]["ph_low"] = None
        data["chemistry"]["ph_high"] = None
    _set_if_present(data["chemistry"], "solvent", overrides.get("solvent"))
    _set_if_present(data["seeding"], "method", overrides.get("seeding_method"))
    _set_if_present(data["refinement"], "censo_enabled", overrides.get("censo_enabled"))
    return RunConfig.model_validate(data)


def write_resolved_config(config: RunConfig, output_dir: Path | None = None) -> Path:
    target_dir = output_dir or config.output_dir
    target_dir.mkdir(parents=True, exist_ok=True)
    path = target_dir / "resolved_config.yaml"
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(_to_yamlable(config.model_dump(mode="python")), handle, sort_keys=False)
    return path


def _set_if_present(data: dict[str, Any], key: str, value: Any) -> None:
    if value is not None:
        data[key] = value


def _to_yamlable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _to_yamlable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_to_yamlable(item) for item in value]
    return value
