from __future__ import annotations

import json
import shutil
from pathlib import Path

from dsvr.config import RunConfig
from dsvr.models import CrestConformerRecord, ThermoRecord, make_thermo_id
from dsvr.parsing.xtb_outputs import parse_xtb_energy, parse_xtb_thermo
from dsvr.runners.subprocess_utils import ExternalToolError, run_command
from dsvr.workflow.provenance import write_jsonl


class XtbUnavailableError(RuntimeError):
    """Raised when xTB is required but not available."""


class XtbExecutionError(RuntimeError):
    """Raised when xTB execution fails unexpectedly."""


def run_xtb_opt(conformer_record: CrestConformerRecord, config: RunConfig) -> Path:
    workdir = xtb_workdir(conformer_record, config)
    workdir.mkdir(parents=True, exist_ok=True)
    input_xyz = _input_xyz_for_conformer(conformer_record)
    target_xyz = workdir / "input.xyz"
    target_xyz.write_text(input_xyz.read_text(encoding="utf-8"), encoding="utf-8")
    command = _base_xtb_command(config, target_xyz, conformer_record)
    command.append("--opt")
    _run_xtb(command, workdir, config, "xtb_opt")
    optimized = workdir / "xtbopt.xyz"
    return optimized if optimized.exists() else target_xyz


def run_xtb_hessian(conformer_record: CrestConformerRecord, config: RunConfig) -> Path:
    workdir = xtb_workdir(conformer_record, config)
    workdir.mkdir(parents=True, exist_ok=True)
    input_xyz = workdir / "xtbopt.xyz"
    if not input_xyz.exists():
        source = _input_xyz_for_conformer(conformer_record)
        input_xyz = workdir / "input.xyz"
        input_xyz.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    command = _base_xtb_command(config, input_xyz, conformer_record)
    command.append("--hess")
    _run_xtb(command, workdir, config, "xtb_hessian")
    return workdir / "hessian"


def run_xtb_thermo(conformer_record: CrestConformerRecord, config: RunConfig) -> ThermoRecord:
    workdir = xtb_workdir(conformer_record, config)
    workdir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []
    try:
        if config.thermo.enabled and config.thermo.xtb_hessian:
            run_xtb_hessian(conformer_record, config)
        if config.thermo.enabled and config.thermo.xtb_thermo:
            command = [
                _xtb_executable(config),
                "thermo",
                "--temp",
                str(config.chemistry.temperature_kelvin),
            ]
            _run_xtb(command, workdir, config, "xtb_thermo")
    except (ExternalToolError, XtbUnavailableError) as exc:
        warnings.append(f"xTB thermo failed: {exc}")

    log_path = _first_existing(
        [
            workdir / "xtb_thermo.out",
            workdir / "thermo.out",
            workdir / "xtb.out",
            workdir / "logs" / "xtb_thermo.log",
        ]
    )
    thermo = parse_xtb_thermo(log_path) if log_path is not None else parse_xtb_energy("")
    free_energy = getattr(thermo, "gibbs_free_energy_kcal_mol", None)
    entropy = getattr(thermo, "entropy_cal_mol_k", None)
    metadata = {
        "xtb": {
            "workdir": str(workdir),
            "logfile": str(log_path) if log_path else None,
            "raw_values": getattr(thermo, "raw_values", None),
            "electronic_energy_kcal_mol": getattr(
                thermo,
                "electronic_energy_kcal_mol",
                None,
            ),
            "enthalpy_kcal_mol": getattr(thermo, "enthalpy_kcal_mol", None),
        }
    }
    record = ThermoRecord(
        id=make_thermo_id(
            conformer_record.id,
            conformer_record.canonical_smiles,
            conformer_record.isomeric_smiles,
            metadata,
        ),
        parent_id=conformer_record.id,
        input_molecule_id=conformer_record.input_molecule_id,
        molname=conformer_record.molname,
        canonical_smiles=conformer_record.canonical_smiles,
        isomeric_smiles=conformer_record.isomeric_smiles,
        molecular_formula=conformer_record.molecular_formula,
        formal_charge=conformer_record.formal_charge,
        explicit_proton_count=conformer_record.explicit_proton_count,
        source_software="xtb",
        source_python_function="dsvr.runners.xtb_runner.run_xtb_thermo",
        output_paths=[workdir / "xtb_thermo.json", workdir / "xtb_provenance.jsonl"],
        warnings=warnings,
        metadata=metadata,
        temperature_kelvin=config.chemistry.temperature_kelvin,
        free_energy_kcal_mol=free_energy,
        entropy_cal_mol_k=entropy,
    )
    _write_thermo_outputs(workdir, record)
    return record


def xtb_workdir(conformer_record: CrestConformerRecord, config: RunConfig) -> Path:
    return (
        config.output_dir
        / "xtb"
        / conformer_record.input_molecule_id
        / (conformer_record.parent_id or "unknown_seed")
        / conformer_record.id
    )


def _base_xtb_command(
    config: RunConfig,
    input_xyz: Path,
    conformer_record: CrestConformerRecord,
) -> list[str]:
    command = [
        _xtb_executable(config),
        str(input_xyz),
        "--gfn",
        str(config.crest.gfn),
        "--chrg",
        str(conformer_record.formal_charge or 0),
        "-P",
        str(config.crest.nproc),
    ]
    if config.chemistry.solvent_model != "none":
        command.extend([f"--{config.chemistry.solvent_model}", config.chemistry.solvent])
    return command


def _run_xtb(command: list[str], workdir: Path, config: RunConfig, name: str) -> None:
    completed = run_command(
        command,
        cwd=workdir,
        timeout_s=None
        if config.crest.walltime_minutes is None
        else config.crest.walltime_minutes * 60,
        log_dir=workdir / "logs",
        command_name=name,
        check=True,
        show_progress=config.logging.tail_subprocess_logs,
    )
    output_path = workdir / f"{name}.out"
    if output_path.exists():
        output_path = workdir / f"{name}.stdout.log"
    output_path.write_text(completed.stdout, encoding="utf-8")


def _xtb_executable(config: RunConfig) -> str:
    path = shutil.which(config.crest.xtb_executable)
    if path is None:
        raise XtbUnavailableError(
            f"xTB executable '{config.crest.xtb_executable}' was not found on PATH."
        )
    return path


def _input_xyz_for_conformer(conformer_record: CrestConformerRecord) -> Path:
    metadata_path = conformer_record.metadata.get("crest", {}).get("workdir")
    if metadata_path is not None:
        candidate = Path(metadata_path) / f"crest_conformer_{conformer_record.crest_index:04d}.xyz"
        if candidate.exists():
            return candidate
    for output_path in conformer_record.output_paths:
        path = Path(output_path)
        if path.name.endswith(".xyz") and path.exists():
            return path
    raise XtbExecutionError(f"No XYZ input found for conformer {conformer_record.id}")


def _write_thermo_outputs(workdir: Path, record: ThermoRecord) -> None:
    (workdir / "xtb_thermo.json").write_text(
        json.dumps(record.model_dump(mode="json"), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_jsonl(workdir / "xtb_provenance.jsonl", [record])


def _first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None
