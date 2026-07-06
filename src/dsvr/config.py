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
Auto3dModel = Literal["AIMNET", "AIMNet2", "ANI2x", "ANI2xt", "auto"]
WorkflowMode = Literal["ligprep_like", "physics_validation", "exhaustive_debug"]
WorkflowProtocol = Literal["default", "auto3d_entropy"]
PopulationScope = Literal["same_formula", "same_charge", "all_approximate"]
EnergyUnit = Literal["kcal/mol"]
LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR"]
QmBackend = Literal["psi4", "pyscf", "none"]
FilteringMode = Literal["conservative", "balanced", "aggressive", "exhaustive"]
CleanupPolicy = Literal["compact", "keep_selected", "debug_all"]
TautomerStrategy = Literal["safe", "normal", "exhaustive"]
TautomerTimeoutFallback = Literal["keep_input_and_canonical"]
OptionalValidationSelection = Literal["top_auto3d_energy"]

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
    max_tautomers_per_protomer: int = 32
    max_tautomer_transforms: int = 256
    tautomer_timeout_seconds: int = 30
    tautomer_strategy: TautomerStrategy = "safe"
    tautomer_remove_bond_stereo: bool = True
    tautomer_remove_sp3_stereo: bool = True
    tautomer_reassign_stereo: bool = True
    max_stereoisomers_per_tautomer: int = 64
    stereo_try_embedding: bool = True
    stereo_only_unassigned: bool = True
    stereo_unique: bool = True
    stereo_random_seed: int = 61453
    fail_on_enumeration_cap: bool = False

    @field_validator(
        "max_protomers_per_molecule",
        "max_tautomers_per_protomer",
        "max_tautomer_transforms",
        "tautomer_timeout_seconds",
        "max_stereoisomers_per_tautomer",
    )
    @classmethod
    def positive_cap(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("enumeration caps must be positive")
        return value


class ProtonationConfig(StrictModel):
    enabled: bool = True
    tool: str = "molscrub"
    mode: str = "plausible"
    max_protomers_per_molecule: int = 4
    keep_input_state: bool = True
    keep_best_per_charge: bool = True
    skip_gen3d_in_molscrub: bool = True
    timeout_seconds_per_molecule: int = 60

    @field_validator("max_protomers_per_molecule", "timeout_seconds_per_molecule")
    @classmethod
    def positive_protonation_limit(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("protonation limits must be positive")
        return value


class TautomerFilteringConfig(StrictModel):
    enabled: bool = True
    tool: str = "auto3d"
    tauto_engine: str = "rdkit"
    optimizing_engine: Auto3dModel = "ANI2xt"
    fallback_optimizing_engine: Auto3dModel = "AIMNET"
    use_gpu: bool = True
    max_rdkit_tautomers_before_auto3d: int = 64
    rdkit_tautomer_timeout_seconds: int = 30
    max_rdkit_transforms: int = 256
    fallback_if_timeout: TautomerTimeoutFallback = "keep_input_and_canonical"
    auto3d_max_confs_per_tautomer: int = 3
    auto3d_patience: int = 100
    tauto_k: int = 3
    tauto_window_kcal_mol: float = 5.0
    keep_input_tautomer: bool = True
    timeout_seconds_per_protomer: int = 600

    @field_validator(
        "max_rdkit_tautomers_before_auto3d",
        "rdkit_tautomer_timeout_seconds",
        "max_rdkit_transforms",
        "auto3d_max_confs_per_tautomer",
        "auto3d_patience",
        "tauto_k",
        "timeout_seconds_per_protomer",
    )
    @classmethod
    def positive_tautomer_filtering_limit(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("tautomer filtering limits must be positive")
        return value

    @field_validator("tauto_window_kcal_mol")
    @classmethod
    def nonnegative_tautomer_window(cls, value: float) -> float:
        if value < 0:
            raise ValueError("tautomer filtering energy window must be non-negative")
        return value


class StereoisomerFilteringConfig(StrictModel):
    enabled: bool = True
    enumerator: str = "rdkit"
    max_stereoisomers_per_tautomer: int = 16
    timeout_seconds_per_tautomer: int = 300
    try_embedding: bool = True
    only_unassigned: bool = True
    collapse_enantiomers_in_achiral_solvent: bool = True
    run_energy_for_enantiomer_representatives_only: bool = True
    stereo_energy_window_kcal_mol: float = 7.0
    keep_top_n_diastereomers: int = 8

    @field_validator(
        "max_stereoisomers_per_tautomer",
        "timeout_seconds_per_tautomer",
        "keep_top_n_diastereomers",
    )
    @classmethod
    def positive_stereoisomer_filtering_limit(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("stereoisomer filtering limits must be positive")
        return value

    @field_validator("stereo_energy_window_kcal_mol")
    @classmethod
    def nonnegative_stereo_window(cls, value: float) -> float:
        if value < 0:
            raise ValueError("stereoisomer filtering energy window must be non-negative")
        return value


class Final3dConfig(StrictModel):
    tool: str = "auto3d"
    k: int = 1
    max_confs: int = 10
    patience: int = 100
    one_conformer_per_variant: bool = True
    optimizing_engine: Auto3dModel = "AIMNET"
    fallback_optimizing_engine: Auto3dModel = "ANI2xt"
    use_gpu: bool = True
    energy_window_kcal_mol: float = 7.0
    timeout_seconds_per_batch: int = 1800

    @field_validator("k", "max_confs", "patience", "timeout_seconds_per_batch")
    @classmethod
    def positive_final_3d_limit(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("final_3d limits must be positive")
        return value

    @field_validator("energy_window_kcal_mol")
    @classmethod
    def nonnegative_final_3d_window(cls, value: float) -> float:
        if value < 0:
            raise ValueError("final_3d energy window must be non-negative")
        return value


class OptionalValidationConfig(StrictModel):
    crest_xtb_enabled: bool = False
    max_variants_per_molecule: int = 5
    selection: OptionalValidationSelection = "top_auto3d_energy"
    xtb_thermo_enabled: bool = False
    crest_entropy_enabled: bool = False
    censo_enabled: bool = False
    keep_raw_xyz: bool = False
    cleanup_policy: CleanupPolicy = "compact"

    @field_validator("max_variants_per_molecule")
    @classmethod
    def positive_optional_validation_limit(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("optional_validation max_variants_per_molecule must be positive")
        return value


class AgentConfig(StrictModel):
    enabled: bool = False
    backend: str = "ollama_codex_cli"
    command: str = "codex --oss -m qwen3.6:35b"
    max_context_chars: int = 12000
    command_timeout_seconds: int = 120
    allowed_tasks: list[str] = Field(
        default_factory=lambda: [
            "classify_failure",
            "suggest_retry_from_menu",
            "summarize_logs",
            "suggest_config_tweak",
        ]
    )
    require_user_approval_for: list[str] = Field(
        default_factory=lambda: [
            "code_patch",
            "science_threshold_change",
            "delete_outputs",
            "rerun_large_job",
        ]
    )
    max_attempts_per_failure: int = 1

    @field_validator("max_context_chars", "command_timeout_seconds", "max_attempts_per_failure")
    @classmethod
    def positive_agent_limit(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("agent limits must be positive")
        return value

    @field_validator("command")
    @classmethod
    def nonempty_agent_command(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("agent command must not be empty")
        return value


class SeedingConfig(StrictModel):
    method: SeederMethod = "etkdg"
    rdkit_num_conformers: int = 30
    rdkit_prune_rms_thresh: float = 0.5
    rdkit_forcefield: RdkitForcefield = "uff"
    auto3d_k: int = 3
    auto3d_model: Auto3dModel = "AIMNet2"
    auto3d_internal_tautomer_stereo_enum: bool = False
    auto3d_mpi_np: int = 4
    auto3d_cpu_workers: int | None = None
    auto3d_memory_gb: int | None = None
    auto3d_capacity: int | None = None
    auto3d_max_confs: int | None = None
    auto3d_patience: int | None = None
    auto3d_threshold: float | None = None
    auto3d_opt_steps: int | None = None
    auto3d_use_gpu: bool = True
    auto3d_allow_rdkit_fallback: bool = True

    @field_validator("rdkit_num_conformers", "auto3d_k", "auto3d_mpi_np")
    @classmethod
    def positive_count(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("conformer counts must be positive")
        return value

    @field_validator(
        "auto3d_cpu_workers",
        "auto3d_memory_gb",
        "auto3d_capacity",
        "auto3d_max_confs",
        "auto3d_patience",
        "auto3d_opt_steps",
    )
    @classmethod
    def optional_positive_auto3d_int(cls, value: int | None) -> int | None:
        if value is not None and value <= 0:
            raise ValueError("Auto3D integer settings must be positive when set")
        return value

    @field_validator("rdkit_prune_rms_thresh")
    @classmethod
    def nonnegative_rms(cls, value: float) -> float:
        if value < 0:
            raise ValueError("rdkit_prune_rms_thresh must be non-negative")
        return value

    @field_validator("auto3d_threshold")
    @classmethod
    def optional_nonnegative_auto3d_threshold(cls, value: float | None) -> float | None:
        if value is not None and value < 0:
            raise ValueError("auto3d_threshold must be non-negative when set")
        return value


class CrestConfig(StrictModel):
    enabled: bool = True
    executable: str = "crest"
    xtb_executable: str = "xtb"
    gfn: int = 2
    ewin_kcal_mol: float = 6.0
    nproc: int = 4
    max_jobs_per_molecule: int = 8
    max_conformers_to_parse: int = 20
    max_conformers_to_keep: int = 5
    keep_raw_xyz: bool = False
    compress_raw_outputs: bool = True
    delete_intermediate_xyz: bool = True
    cleanup_patterns: list[str] = Field(
        default_factory=lambda: ["*.tmp", "coord.*", "struc*.xyz", "trial*.xyz"]
    )
    walltime_minutes: int | None = 30
    command_template: str | None = None
    extra_args: list[str] = Field(default_factory=list)

    @field_validator(
        "gfn",
        "nproc",
        "max_jobs_per_molecule",
        "max_conformers_to_parse",
        "max_conformers_to_keep",
    )
    @classmethod
    def positive_int(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("CREST positive integer settings must be positive")
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


class XtbPrefilterConfig(StrictModel):
    enabled: bool = False
    optimize: bool = True
    hessian: bool = False
    solvent_model: SolventModel = "alpb"
    solvent: str = "water"
    gfn: int = 2
    max_variants_per_molecule: int = 16
    keep_within_kcal_mol: float = 10.0
    keep_top_n_per_molecule: int = 8
    keep_top_n_per_charge: int = 3
    keep_top_n_per_formula: int = 3
    timeout_seconds_per_variant: int = 300
    nproc: int = 4

    @field_validator(
        "gfn",
        "max_variants_per_molecule",
        "keep_top_n_per_molecule",
        "keep_top_n_per_charge",
        "keep_top_n_per_formula",
        "timeout_seconds_per_variant",
        "nproc",
    )
    @classmethod
    def positive_prefilter_int(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("xTB prefilter integer limits must be positive")
        return value

    @field_validator("keep_within_kcal_mol")
    @classmethod
    def nonnegative_prefilter_window(cls, value: float) -> float:
        if value < 0:
            raise ValueError("xTB prefilter energy window must be non-negative")
        return value


class ThermoConfig(StrictModel):
    enabled: bool = True
    xtb_hessian: bool = True
    xtb_thermo: bool = True
    rrho_cutoff: float = 100.0
    population_scope: PopulationScope = "same_formula"
    energy_unit: EnergyUnit = "kcal/mol"
    max_variants_per_molecule: int = 5
    max_conformers_per_variant: int = 3

    @field_validator("rrho_cutoff")
    @classmethod
    def nonnegative_rrho_cutoff(cls, value: float) -> float:
        if value < 0:
            raise ValueError("rrho_cutoff must be non-negative")
        return value

    @field_validator("max_variants_per_molecule", "max_conformers_per_variant")
    @classmethod
    def positive_thermo_limit(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("thermo candidate limits must be positive")
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
    max_candidates_for_refinement: int = 3

    @field_validator("max_candidates_for_refinement")
    @classmethod
    def positive_refinement_count(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("max_candidates_for_refinement must be positive")
        return value


class QmConfig(StrictModel):
    enabled: bool = False
    backend: QmBackend = "none"
    max_candidates: int = 3

    @field_validator("max_candidates")
    @classmethod
    def positive_qm_candidates(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("qm.max_candidates must be positive")
        return value


class LoggingConfig(StrictModel):
    level: LogLevel = "INFO"
    save_subprocess_logs: bool = True
    tail_subprocess_logs: bool = True


class VariantFilteringConfig(StrictModel):
    enabled: bool = True
    mode: FilteringMode = "balanced"
    max_variants_before_3d_per_molecule: int = 64
    max_variants_after_cheap_score_per_molecule: int = 24
    max_variants_for_xtb_prefilter_per_molecule: int = 16
    max_variants_for_crest_per_molecule: int = 8
    max_seeds_per_variant: int = 2
    keep_original_state: bool = True
    keep_best_per_charge_state: bool = True
    keep_best_per_formula: bool = True
    keep_best_per_protomer: bool = True
    keep_best_per_tautomer_family: bool = True
    collapse_enantiomers: bool = True
    absolute_penalty_cutoff: float = 12.0
    relative_penalty_cutoff: float = 7.0
    min_variants_to_keep: int = 3
    rescue_rules_enabled: bool = True

    @field_validator(
        "max_variants_before_3d_per_molecule",
        "max_variants_after_cheap_score_per_molecule",
        "max_variants_for_xtb_prefilter_per_molecule",
        "max_variants_for_crest_per_molecule",
        "max_seeds_per_variant",
        "min_variants_to_keep",
    )
    @classmethod
    def positive_budget(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("variant filtering budgets must be positive")
        return value

    @field_validator("absolute_penalty_cutoff", "relative_penalty_cutoff")
    @classmethod
    def nonnegative_cutoff(cls, value: float) -> float:
        if value < 0:
            raise ValueError("variant filtering cutoffs must be non-negative")
        return value


class StereoFilteringConfig(StrictModel):
    collapse_enantiomers_in_achiral_solvent: bool = True
    solvent_is_chiral: bool = False
    run_crest_for_enantiomer_pairs_once: bool = True
    keep_mapping_to_all_stereo_outputs: bool = True


class TimeoutConfig(StrictModel):
    protomer_seconds_per_molecule: int = 60
    tautomer_seconds_per_protomer: int = 30
    stereo_seconds_per_tautomer: int = 30
    etkdg_seconds_per_variant: int = 60
    auto3d_seconds_per_batch: int = 600
    xtb_prefilter_seconds_per_variant: int = 300
    crest_seconds_per_variant: int = 1800
    thermo_seconds_per_conformer: int = 600

    @field_validator("*")
    @classmethod
    def positive_timeout(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("timeouts must be positive")
        return value


class DiskConfig(StrictModel):
    cleanup_policy: CleanupPolicy = "compact"
    keep_raw_xyz: bool = False
    compress_raw_outputs: bool = True
    delete_intermediate_xyz: bool = True
    max_run_dir_gb: float = 20.0
    max_molecule_dir_gb: float = 3.0
    max_xyz_files_per_molecule: int = 500
    fail_on_disk_limit: bool = True

    @field_validator("max_run_dir_gb", "max_molecule_dir_gb")
    @classmethod
    def positive_disk_limit(cls, value: float) -> float:
        if value <= 0:
            raise ValueError("disk GB limits must be positive")
        return value

    @field_validator("max_xyz_files_per_molecule")
    @classmethod
    def positive_xyz_limit(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("max_xyz_files_per_molecule must be positive")
        return value


class ErrorHandlingConfig(StrictModel):
    fail_fast: bool = False
    keep_fallback_parent_state: bool = True
    skip_failed_molecule: bool = True
    skip_failed_variant: bool = True
    retry_auto3d_cpu_on_gpu_failure: bool = True
    retry_auto3d_smaller_batch_on_batch_failure: bool = True
    reduce_tautomer_cap_on_timeout: bool = True
    reduce_stereo_cap_on_timeout: bool = True


class RunConfig(StrictModel):
    description: str | None = None
    workflow_mode: WorkflowMode = "ligprep_like"
    protocol: WorkflowProtocol = "default"
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
    protonation: ProtonationConfig = Field(default_factory=ProtonationConfig)
    tautomer_filtering: TautomerFilteringConfig = Field(default_factory=TautomerFilteringConfig)
    stereoisomer_filtering: StereoisomerFilteringConfig = Field(default_factory=StereoisomerFilteringConfig)
    final_3d: Final3dConfig = Field(default_factory=Final3dConfig)
    optional_validation: OptionalValidationConfig = Field(default_factory=OptionalValidationConfig)
    agent: AgentConfig = Field(default_factory=AgentConfig)
    seeding: SeedingConfig = Field(default_factory=SeedingConfig)
    crest: CrestConfig = Field(default_factory=CrestConfig)
    xtb_prefilter: XtbPrefilterConfig = Field(default_factory=XtbPrefilterConfig)
    thermo: ThermoConfig = Field(default_factory=ThermoConfig)
    refinement: OptionalRefinementConfig = Field(default_factory=OptionalRefinementConfig)
    qm: QmConfig = Field(default_factory=QmConfig)
    logging: LoggingConfig = Field(default_factory=LoggingConfig)
    variant_filtering: VariantFilteringConfig = Field(default_factory=VariantFilteringConfig)
    stereo_filtering: StereoFilteringConfig = Field(default_factory=StereoFilteringConfig)
    timeouts: TimeoutConfig = Field(default_factory=TimeoutConfig)
    disk: DiskConfig = Field(default_factory=DiskConfig)
    error_handling: ErrorHandlingConfig = Field(default_factory=ErrorHandlingConfig)

    @field_validator("max_workers")
    @classmethod
    def positive_max_workers(cls, value: int) -> int:
        if value <= 0:
            raise ValueError("max_workers must be positive")
        return value

    @model_validator(mode="before")
    @classmethod
    def migrate_legacy_workflow_mode(cls, data: Any) -> Any:
        if not isinstance(data, dict):
            return data
        chemistry = data.get("chemistry")
        if isinstance(chemistry, dict) and "workflow_mode" in chemistry:
            data = dict(data)
            chemistry = dict(chemistry)
            legacy_mode = chemistry.pop("workflow_mode")
            data["chemistry"] = chemistry
            if "workflow_mode" not in data:
                data["workflow_mode"] = legacy_mode
        return data

    @model_validator(mode="after")
    def validate_ligprep_like_settings(self) -> RunConfig:
        if self.workflow_mode == "ligprep_like":
            if not self.final_3d.one_conformer_per_variant:
                raise ValueError(
                    "ligprep_like mode requires final_3d.one_conformer_per_variant=true"
                )
            if self.tautomer_filtering.timeout_seconds_per_protomer <= 0:
                raise ValueError("tautomer enumeration timeout must be finite and positive")
            if self.stereoisomer_filtering.timeout_seconds_per_tautomer <= 0:
                raise ValueError("stereoisomer enumeration timeout must be finite and positive")
        return self


def load_config(path: Path) -> RunConfig:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if path.name == "resolved_config.yaml":
        output_dir = data.get("output_dir")
        if output_dir is not None and not Path(output_dir).is_absolute():
            data["output_dir"] = path.parent
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
    _set_if_present(
        data["optional_validation"],
        "crest_xtb_enabled",
        overrides.get("crest_xtb_enabled"),
    )
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
