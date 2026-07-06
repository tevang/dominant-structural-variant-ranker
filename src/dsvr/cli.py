from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Any

import typer
from rich.console import Console
from rich.table import Table

from dsvr.agent.action_menu import deterministic_action_for_failure
from dsvr.agent.bug_package import build_bug_package
from dsvr.agent.local_qwen import LocalAgentResult, run_local_diagnostic_agent
from dsvr.chemistry.conformers_auto3d import generate_auto3d_seeds
from dsvr.chemistry.conformers_rdkit import generate_rdkit_seeds, read_stereo_sdf
from dsvr.chemistry.protonation import generate_protomer_candidates
from dsvr.chemistry.stereochemistry import enumerate_stereoisomers, read_tautomers_sdf
from dsvr.chemistry.tautomers import enumerate_tautomers, read_protomers_sdf
from dsvr.config import RunConfig, SeederMethod, load_config, merge_cli_overrides
from dsvr.io.read_inputs import InputFormat, read_molecules, validate_input_file
from dsvr.io.write_outputs import write_final_ranked_outputs, write_json
from dsvr.models import CrestConformerRecord, RankedVariantRecord, ThermoRecord
from dsvr.ranking.population import (
    compute_delta_g_and_populations,
    load_rankable_records,
    write_ranked_outputs,
)
from dsvr.reporting.markdown import write_run_report
from dsvr.runners.auto3d_runner import Auto3DExecutionError, Auto3DUnavailableError
from dsvr.runners.crest_runner import (
    CrestExecutionError,
    CrestUnavailableError,
    read_seed_sdf,
    run_crest_for_seed,
)
from dsvr.runners.molscrub_runner import MolscrubUnavailableError
from dsvr.runners.psi4_runner import (
    Psi4UnavailableError,
    load_ranked_for_qm,
    rescore_top_ranked_with_psi4,
)
from dsvr.runners.pyscf_runner import PySCFUnavailableError, rescore_top_ranked_with_pyscf
from dsvr.runners.xtb_runner import XtbExecutionError, XtbUnavailableError, run_xtb_thermo
from dsvr.utils.tool_check import check_tools
from dsvr.workflow.engine import run_workflow
from dsvr.workflow.status import run_status
from dsvr.workflow.steps import planned_steps

app = typer.Typer(help="Rank pH- and solvent-dependent structural variants of small molecules.")
agent_app = typer.Typer(help="Experimental opt-in local diagnostic agent commands.")
app.add_typer(agent_app, name="agent")
console = Console()
LOG_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR"}


def _resolve_log_level(*, verbose: bool, quiet: bool, log_level: str | None) -> str | None:
    if log_level is not None:
        normalized = log_level.upper()
        if normalized not in LOG_LEVELS:
            raise typer.BadParameter("log-level must be one of DEBUG, INFO, WARNING, ERROR")
        return normalized
    if verbose:
        return "DEBUG"
    if quiet:
        return "ERROR"
    return None


def _global_options(ctx: typer.Context) -> dict[str, Any]:
    return dict(ctx.obj or {})


def _apply_runtime_options(
    config: RunConfig,
    *,
    globals_: dict[str, Any],
    dry_run: bool | None = None,
    overwrite: bool | None = None,
    resume: bool | None = None,
) -> RunConfig:
    data = config.model_dump(mode="python")
    if dry_run is not None:
        data["dry_run"] = dry_run
    elif globals_.get("dry_run"):
        data["dry_run"] = True

    overwrite_value = overwrite if overwrite is not None else globals_.get("overwrite")
    if overwrite_value is not None:
        data["overwrite"] = overwrite_value

    resume_value = resume if resume is not None else globals_.get("resume")
    if resume_value is not None:
        data["resume"] = resume_value

    if globals_.get("log_level") is not None:
        data["logging"]["level"] = globals_["log_level"]
    return RunConfig.model_validate(data)


def _validation_report_path(out: Path) -> Path:
    if out.suffix.lower() == ".json":
        return out
    return out / "validation_report.json"


def _load_agent_enabled_run_config(run_dir: Path) -> RunConfig:
    resolved_config = run_dir / "resolved_config.yaml"
    if not resolved_config.exists():
        raise typer.BadParameter(
            f"{run_dir} does not contain resolved_config.yaml. Cannot run agent diagnostics."
        )
    config = load_config(resolved_config)
    data = config.model_dump(mode="python")
    data["output_dir"] = run_dir
    config = RunConfig.model_validate(data)
    if not config.agent.enabled:
        raise typer.BadParameter(
            "Local diagnostic agent is disabled. Set agent.enabled=true in resolved_config.yaml "
            "or rerun with an enabled config to use dsvr agent commands."
        )
    return config


def _run_agent_diagnostic(
    run_dir: Path,
    *,
    latest: bool,
    task_key: str,
    task: str,
) -> dict[str, Any]:
    config = _load_agent_enabled_run_config(run_dir)
    if task_key not in config.agent.allowed_tasks:
        raise typer.BadParameter(f"Agent task is not allowed by config: {task_key}")
    package = build_bug_package(run_dir, config, latest=latest)
    deterministic_action = deterministic_action_for_failure(package.failure.get("failure_kind"))
    result = run_local_diagnostic_agent(
        agent=config.agent,
        task=task,
        bug_context=package.prompt_context,
    )
    return _agent_payload(
        run_dir=run_dir,
        package_path=package.path,
        failure=package.failure,
        deterministic_action=deterministic_action,
        result=result,
    )


def _agent_payload(
    *,
    run_dir: Path,
    package_path: Path,
    failure: dict[str, Any],
    deterministic_action: str,
    result: LocalAgentResult,
) -> dict[str, Any]:
    decision = result.decision
    return {
        "run_dir": str(run_dir),
        "bug_package": str(package_path),
        "failure_kind": failure.get("failure_kind", "UNKNOWN"),
        "stage": failure.get("stage", ""),
        "item_id": failure.get("item_id"),
        "deterministic_action": deterministic_action,
        "agent_available": result.available,
        "agent_action": decision.action,
        "agent_output_valid": decision.valid,
        "reasons": decision.reasons,
        "config_tweak": decision.config_tweak,
        "error": result.error,
    }


def _print_json_payload(payload: dict[str, Any]) -> None:
    typer.echo(json.dumps(payload, indent=2, sort_keys=True, default=str))


@app.callback()
def main(
    ctx: typer.Context,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", help="Enable DEBUG-level console/log output."),
    ] = False,
    quiet: Annotated[
        bool,
        typer.Option("--quiet", help="Reduce console/log output to ERROR-level messages."),
    ] = False,
    log_level: Annotated[
        str | None,
        typer.Option(
            "--log-level",
            help="Explicit log level: DEBUG, INFO, WARNING, or ERROR.",
        ),
    ] = None,
    config_path: Annotated[
        Path | None,
        typer.Option(
            "--config",
            "-c",
            exists=True,
            dir_okay=False,
            help="Global YAML workflow config. `dsvr run --config` overrides this.",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Global dry-run flag for commands that support it."),
    ] = False,
    overwrite: Annotated[
        bool | None,
        typer.Option(
            "--overwrite/--no-overwrite",
            help="Global overwrite policy for commands that support it.",
        ),
    ] = None,
    resume: Annotated[
        bool | None,
        typer.Option(
            "--resume/--no-resume",
            help="Global resume policy for commands that support it.",
        ),
    ] = None,
) -> None:
    """Dominant Structural Variant Ranker CLI."""
    if verbose and quiet:
        raise typer.BadParameter("--verbose and --quiet are mutually exclusive")
    resolved_log_level = _resolve_log_level(
        verbose=verbose,
        quiet=quiet,
        log_level=log_level,
    )
    ctx.obj = {
        "config_path": config_path,
        "dry_run": dry_run,
        "overwrite": overwrite,
        "resume": resume,
        "log_level": resolved_log_level,
    }


@app.command()
def version() -> None:
    """Print package version."""
    from dsvr import __version__

    console.print(__version__)


@app.command()
def doctor(
    output_dir: Annotated[
        Path,
        typer.Option("--output-dir", help="Directory used for writability and disk checks."),
    ] = Path("runs/dsvr"),
    json_report: Annotated[
        bool,
        typer.Option("--json", help="Write a machine-readable JSON doctor report."),
    ] = False,
    json_out: Annotated[
        Path,
        typer.Option("--json-out", help="Path for --json doctor report."),
    ] = Path("doctor_report.json"),
    strict: Annotated[
        bool,
        typer.Option("--strict", help="Exit nonzero if required default-workflow checks fail."),
    ] = False,
) -> None:
    """Check optional external tools without importing them at package import time."""
    table = Table(title="DSVR external tool check")
    table.add_column("Tool")
    table.add_column("Kind")
    table.add_column("Required")
    table.add_column("Status")
    table.add_column("Version")
    table.add_column("Detail")
    statuses = check_tools(output_dir=output_dir)
    required_missing = [
        item.name for item in statuses if item.required and not item.available
    ]
    payload = {
        "ok": not required_missing,
        "strict_failure_count": len(required_missing),
        "required_missing": required_missing,
        "checks": [item.model_dump(mode="json") for item in statuses],
    }
    for item in statuses:
        table.add_row(
            item.name,
            item.kind,
            "yes" if item.required else "optional",
            "ok" if item.available else "missing",
            item.version or "",
            item.detail,
        )
    console.print(table)
    if json_report:
        write_json(json_out, payload)
        console.print(f"Wrote doctor JSON report to [bold]{json_out}[/bold]")
    if strict and payload["required_missing"]:
        missing = ", ".join(payload["required_missing"])
        console.print(f"[red]Required checks failed: {missing}[/red]")
        raise typer.Exit(1)


@app.command()
def inspect(
    input_path: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
) -> None:
    """Inspect input molecules without running external tools."""
    molecules = read_molecules(input_path)
    table = Table(title=f"Input molecules: {input_path}")
    table.add_column("Index")
    table.add_column("Name")
    table.add_column("Format")
    table.add_column("SMILES")
    for mol in molecules:
        table.add_row(str(mol.index), mol.name, mol.input_format, mol.smiles or "")
    console.print(table)


@agent_app.command(name="diagnose")
def agent_diagnose(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
    latest: Annotated[
        bool,
        typer.Option("--latest/--no-latest", help="Diagnose the latest recorded failure."),
    ] = True,
) -> None:
    """Build a compact bug package and run bounded local failure diagnostics."""
    payload = _run_agent_diagnostic(
        run_dir,
        latest=latest,
        task_key="classify_failure",
        task="Summarize and classify the latest compact failure package.",
    )
    _print_json_payload(payload)


@agent_app.command(name="suggest-retry")
def agent_suggest_retry(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
    latest: Annotated[
        bool,
        typer.Option("--latest/--no-latest", help="Suggest action for the latest recorded failure."),
    ] = True,
) -> None:
    """Ask the local agent to choose one action from the fixed retry menu."""
    payload = _run_agent_diagnostic(
        run_dir,
        latest=latest,
        task_key="suggest_retry_from_menu",
        task="Choose exactly one retry or skip action from the allowed action menu.",
    )
    _print_json_payload(payload)


@agent_app.command(name="summarize")
def agent_summarize(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
    latest: Annotated[
        bool,
        typer.Option("--latest/--no-latest", help="Summarize the latest recorded failure."),
    ] = True,
) -> None:
    """Ask the local agent for a constrained summary of the compact bug package."""
    payload = _run_agent_diagnostic(
        run_dir,
        latest=latest,
        task_key="summarize_logs",
        task="Summarize the compact failure package and choose one allowed action.",
    )
    _print_json_payload(payload)


@app.command(name="status")
def status_command(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
) -> None:
    """Show current or last known workflow status for a run directory."""
    status = run_status(run_dir)
    table = Table(title=f"DSVR status: {run_dir}")
    table.add_column("Field")
    table.add_column("Value")
    current = _format_current_molecule(status)
    for key, value in (
        ("last_completed_stage", status.get("last_completed_stage")),
        ("active_stage", status.get("active_stage")),
        ("last_stage", status.get("last_stage")),
        ("last_status", status.get("last_status")),
        ("current_molecule", current),
        ("active_command", status.get("active_command")),
        ("disk_usage_mb", f"{float(status.get('disk_usage_mb', 0.0)):.2f}"),
        ("xyz_file_count", status.get("xyz_file_count")),
        ("resume_possible", status.get("resume_possible")),
    ):
        table.add_row(key, "" if value is None else str(value))
    console.print(table)

    rows = list(status.get("counts_by_stage", []))
    if rows:
        counts = Table(title="Counts by stage")
        columns = [
            "stage",
            "status",
            "generated_count",
            "selected_count",
            "rejected_count",
            "timeout_count",
        ]
        for column in columns:
            counts.add_column(column)
        for row in rows:
            counts.add_row(*(str(row.get(column, "")) for column in columns))
        console.print(counts)
    else:
        console.print("Stage counts:")
        for stage, count in dict(status.get("stage_counts", {})).items():
            console.print(f"- {stage}: {count}")

    _print_diagnostics("Latest warnings", status.get("latest_warnings", []))
    _print_diagnostics("Latest failures", status.get("latest_failures", []))
    _print_diagnostics("Latest recovery failures", status.get("latest_recovery_failures", []))
    if status.get("molecule_state_count"):
        console.print(f"Molecule state files: {status.get('molecule_state_count')}")
    if status.get("latest_log_tail"):
        console.print("Latest log tail:")
        console.print(str(status["latest_log_tail"]))


def _format_current_molecule(status: dict[str, Any]) -> str:
    name = status.get("current_molecule") or ""
    index = status.get("current_molecule_index")
    total = status.get("current_molecule_total")
    if index is not None and total is not None:
        return f"{index}/{total} {name}".strip()
    return str(name)


def _print_diagnostics(title: str, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    console.print(f"{title}:")
    for row in rows:
        stage = row.get("stage", "")
        message = row.get("message", "")
        console.print(f"- {stage}: {message}" if stage else f"- {message}")


@app.command(name="resume")
def resume_command(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
) -> None:
    """Resume a deterministic DSVR workflow run from its run directory."""
    resolved_config = run_dir / "resolved_config.yaml"
    if not resolved_config.exists():
        raise typer.BadParameter(
            f"{run_dir} does not contain resolved_config.yaml. Cannot resume deterministically."
        )
    config = load_config(resolved_config)
    data = config.model_dump(mode="python")
    data["output_dir"] = run_dir
    data["resume"] = True
    data["overwrite"] = False
    data["dry_run"] = False
    config = RunConfig.model_validate(data)
    result = run_workflow(config=config)
    console.print(f"Resumed workflow outputs in [bold]{result.outdir}[/bold]")
    console.print(f"Molecules: {result.molecule_count}")


@app.command(name="validate-input")
def validate_input(
    input_path: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    input_format: Annotated[
        InputFormat,
        typer.Option("--format", help="Input format override."),
    ] = "auto",
    out: Annotated[
        Path,
        typer.Option("--out", help="Validation report JSON path."),
    ] = Path("validation_report.json"),
    deduplicate: Annotated[
        bool,
        typer.Option(
            "--deduplicate/--no-deduplicate",
            help="Deduplicate by canonical isomeric SMILES.",
        ),
    ] = True,
) -> None:
    """Validate input molecules and write a JSON report plus invalid_inputs.csv."""
    report_path = _validation_report_path(out)
    invalid_path = report_path.parent / "invalid_inputs.csv"
    molecules, invalid_records = validate_input_file(
        input_path,
        input_format=input_format,
        deduplicate=deduplicate,
        invalid_output_path=invalid_path,
    )
    report = {
        "input_path": str(input_path),
        "input_format": input_format,
        "deduplicate": deduplicate,
        "valid_count": len(molecules),
        "invalid_count": len(invalid_records),
        "invalid_inputs_csv": str(invalid_path) if invalid_records else None,
        "molecules": [
            {
                "input_id": molecule.input_id,
                "molname": molecule.molname,
                "source_format": molecule.source_format,
                "original_smiles": molecule.original_smiles,
                "canonical_smiles": molecule.canonical_smiles,
                "isomeric_smiles": molecule.isomeric_smiles,
                "warnings": molecule.warnings,
            }
            for molecule in molecules
        ],
        "invalid_records": invalid_records,
    }
    write_json(report_path, report)
    console.print(f"Valid molecules: {len(molecules)}")
    console.print(f"Invalid records: {len(invalid_records)}")
    console.print(f"Wrote validation report to [bold]{report_path}[/bold]")


@app.command(name="enumerate-protomers")
def enumerate_protomers(
    input_path: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    ph: Annotated[
        float,
        typer.Option("--ph", help="Candidate-generation pH for molscrub."),
    ] = 7.0,
    solvent: Annotated[
        str,
        typer.Option("--solvent", help="Solvent label recorded in protomer metadata."),
    ] = "water",
    out: Annotated[
        Path,
        typer.Option("--out", help="Output run directory."),
    ] = Path("runs/test_protomers"),
    input_format: Annotated[
        InputFormat,
        typer.Option("--format", help="Input format override."),
    ] = "auto",
) -> None:
    """Generate pH/protomer candidate states with molscrub."""
    config = RunConfig(
        input_path=input_path,
        input_format=input_format,
        output_dir=out,
        chemistry={"ph": ph, "solvent": solvent},
    )
    molecules = read_molecules(
        input_path,
        input_format=input_format,
        invalid_output_path=out / "invalid_inputs.csv",
    )
    all_records = []
    try:
        for molecule in molecules:
            all_records.extend(generate_protomer_candidates(molecule, config))
    except MolscrubUnavailableError as exc:
        raise typer.BadParameter(str(exc)) from exc

    report = {
        "input_path": str(input_path),
        "candidate_generation_ph": ph,
        "solvent": solvent,
        "protomer_count": len(all_records),
        "scope_warning": (
            "molscrub pH influence is candidate generation only; no rigorous pH "
            "population prediction is implied."
        ),
        "protomer_ids": [record.id for record in all_records],
    }
    write_json(out / "enumeration" / "protomers" / "protomer_report.json", report)
    console.print(f"Generated protomer candidates: {len(all_records)}")
    console.print("pH scope: candidate generation only; not rigorous pH population prediction")
    console.print(f"Wrote outputs to [bold]{out / 'enumeration' / 'protomers'}[/bold]")


@app.command(name="enumerate-tautomers")
def enumerate_tautomers_command(
    protomers_sdf: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    out: Annotated[
        Path,
        typer.Option("--out", help="Output run directory."),
    ] = Path("runs/test_tautomers"),
    max_tautomers: Annotated[
        int | None,
        typer.Option("--max-tautomers", help="Override max tautomers per protomer."),
    ] = None,
    max_transforms: Annotated[
        int | None,
        typer.Option("--max-transforms", help="Override RDKit max tautomer transforms."),
    ] = None,
    timeout_seconds: Annotated[
        int | None,
        typer.Option("--timeout", help="Tautomer enumeration timeout per protomer in seconds."),
    ] = None,
    strategy: Annotated[
        str | None,
        typer.Option("--strategy", help="Tautomer strategy: safe, normal, or exhaustive."),
    ] = None,
) -> None:
    """Enumerate RDKit tautomer candidate states from protomer SDF."""
    config_kwargs = {"input_path": protomers_sdf, "output_dir": out}
    enumeration_config: dict[str, Any] = {}
    if max_tautomers is not None:
        enumeration_config["max_tautomers_per_protomer"] = max_tautomers
    if max_transforms is not None:
        enumeration_config["max_tautomer_transforms"] = max_transforms
    if timeout_seconds is not None:
        enumeration_config["tautomer_timeout_seconds"] = timeout_seconds
    if strategy is not None:
        if strategy not in {"safe", "normal", "exhaustive"}:
            raise typer.BadParameter("strategy must be one of: safe, normal, exhaustive")
        enumeration_config["tautomer_strategy"] = strategy
    if enumeration_config:
        config_kwargs["enumeration"] = enumeration_config
    config = RunConfig(**config_kwargs)
    protomers = read_protomers_sdf(protomers_sdf)
    all_records = []
    for protomer in protomers:
        all_records.extend(enumerate_tautomers(protomer, config))

    report = {
        "input_path": str(protomers_sdf),
        "protomer_count": len(protomers),
        "tautomer_count": len(all_records),
        "scope_warning": (
            "RDKit tautomer enumeration is candidate generation only; no tautomer "
            "stability ranking is implied. CREST/xTB ranking occurs later."
        ),
        "timeout_count": sum(
            1
            for record in all_records
            if any("tautomer enumeration timeout" in warning for warning in record.warnings)
        ),
        "tautomer_ids": [record.id for record in all_records],
    }
    write_json(out / "enumeration" / "tautomers" / "tautomer_report.json", report)
    console.print(f"Enumerated tautomer candidates: {len(all_records)}")
    console.print("Scope: candidate generation only; not tautomer stability ranking")
    console.print(f"Wrote outputs to [bold]{out / 'enumeration' / 'tautomers'}[/bold]")


@app.command(name="enumerate-stereo")
def enumerate_stereo_command(
    tautomers_sdf: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    out: Annotated[
        Path,
        typer.Option("--out", help="Output run directory."),
    ] = Path("runs/test_stereo"),
    max_isomers: Annotated[
        int | None,
        typer.Option("--max-isomers", help="Override max stereoisomers per tautomer."),
    ] = None,
    enumerate_assigned: Annotated[
        bool,
        typer.Option(
            "--enumerate-assigned/--preserve-assigned",
            help="Enumerate already assigned stereocenters instead of preserving them.",
        ),
    ] = False,
) -> None:
    """Enumerate explicit RDKit stereoisomer candidate states from tautomer SDF."""
    enumeration_config = {"stereo_only_unassigned": not enumerate_assigned}
    if max_isomers is not None:
        enumeration_config["max_stereoisomers_per_tautomer"] = max_isomers
    config = RunConfig(
        input_path=tautomers_sdf,
        output_dir=out,
        enumeration=enumeration_config,
    )
    tautomers = read_tautomers_sdf(tautomers_sdf)
    all_records = []
    for tautomer in tautomers:
        all_records.extend(enumerate_stereoisomers(tautomer, config))

    report = {
        "input_path": str(tautomers_sdf),
        "tautomer_count": len(tautomers),
        "stereoisomer_count": len(all_records),
        "try_embedding_scope": (
            "tryEmbedding is a heuristic filter and can be computationally expensive."
        ),
        "preserve_assigned_stereo": not enumerate_assigned,
        "stereo_ids": [record.id for record in all_records],
    }
    write_json(out / "enumeration" / "stereoisomers" / "stereo_report.json", report)
    console.print(f"Enumerated stereoisomer candidates: {len(all_records)}")
    console.print("tryEmbedding is heuristic and can be computationally expensive")
    console.print(f"Wrote outputs to [bold]{out / 'enumeration' / 'stereoisomers'}[/bold]")


@app.command(name="seed-etkdg")
def seed_etkdg_command(
    stereo_sdf: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    out: Annotated[
        Path,
        typer.Option("--out", help="Output run directory."),
    ] = Path("runs/test_seeds"),
    num_conformers: Annotated[
        int,
        typer.Option("--num-conformers", help="Number of RDKit conformers to request."),
    ] = 20,
) -> None:
    """Generate RDKit ETKDG seed conformers and XYZ files for xTB/CREST."""
    config = RunConfig(
        input_path=stereo_sdf,
        output_dir=out,
        seeding={"rdkit_num_conformers": num_conformers},
        disk={"keep_raw_xyz": True},
    )
    stereo_records = read_stereo_sdf(stereo_sdf)
    all_records = []
    for stereo_record in stereo_records:
        all_records.extend(generate_rdkit_seeds(stereo_record, config))
    report = {
        "input_path": str(stereo_sdf),
        "stereo_count": len(stereo_records),
        "seed_count": len(
            [record for record in all_records if record.embedding_status == "success"]
        ),
        "failure_count": len(
            [record for record in all_records if record.embedding_status == "failed"]
        ),
        "seed_ids": [record.id for record in all_records],
    }
    write_json(out / "seeding" / "rdkit" / "seed_report.json", report)
    console.print(f"Generated RDKit seed records: {len(all_records)}")
    console.print(f"Wrote outputs to [bold]{out / 'seeding' / 'rdkit'}[/bold]")


@app.command(name="seed-auto3d")
def seed_auto3d_command(
    stereo_sdf: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    out: Annotated[
        Path,
        typer.Option("--out", help="Output run directory."),
    ] = Path("runs/test_auto3d"),
    k: Annotated[
        int,
        typer.Option("--k", help="Number of Auto3D conformers per input."),
    ] = 5,
    model: Annotated[
        str,
        typer.Option("--model", help="Auto3D model: AIMNet2, ANI2x, ANI2xt, or auto."),
    ] = "AIMNet2",
    internal_enum: Annotated[
        bool,
        typer.Option(
            "--internal-enum/--no-internal-enum",
            help=(
                "Allow Auto3D internal tautomer/stereo enumeration. Disabled by default "
                "to avoid double enumeration."
            ),
        ),
    ] = False,
) -> None:
    """Generate optional Auto3D seed conformers from stereoisomer SDF input."""
    config = RunConfig(
        input_path=stereo_sdf,
        output_dir=out,
        seeding={
            "method": "auto3d",
            "auto3d_k": k,
            "auto3d_model": model,
            "auto3d_internal_tautomer_stereo_enum": internal_enum,
        },
    )
    stereo_records = read_stereo_sdf(stereo_sdf)
    try:
        records = generate_auto3d_seeds(stereo_records, config)
    except (Auto3DUnavailableError, Auto3DExecutionError) as exc:
        raise typer.BadParameter(str(exc)) from exc

    report = {
        "input_path": str(stereo_sdf),
        "stereo_count": len(stereo_records),
        "seed_count": len(records),
        "auto3d_k": config.seeding.auto3d_k,
        "auto3d_model": config.seeding.auto3d_model,
        "internal_tautomer_stereo_enum": (
            config.seeding.auto3d_internal_tautomer_stereo_enum
        ),
        "lineage_scope": (
            "auto3d_internal_enum"
            if config.seeding.auto3d_internal_tautomer_stereo_enum
            else "post_stereo_seed"
        ),
        "double_enumeration_warning": (
            "Auto3D internal tautomer/stereo enumeration is disabled by default to "
            "preserve explicit RDKit protomer-tautomer-stereo lineage."
        ),
        "seed_ids": [record.id for record in records],
    }
    write_json(out / "seeding" / "auto3d" / "auto3d_report.json", report)
    console.print(f"Generated Auto3D seed records: {len(records)}")
    if internal_enum:
        console.print(
            "Lineage scope: Auto3D internal enumeration; exact protomer-tautomer-stereo "
            "lineage may be less controlled"
        )
    else:
        console.print("Auto3D internal tautomer/stereo enumeration disabled")
    console.print(f"Wrote outputs to [bold]{out / 'seeding' / 'auto3d'}[/bold]")


@app.command(name="run-crest")
def run_crest_command(
    seeds_sdf: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    out: Annotated[
        Path,
        typer.Option("--out", help="Output run directory."),
    ] = Path("runs/test_crest"),
    solvent: Annotated[
        str,
        typer.Option("--solvent", help="Solvent passed to CREST/xTB where supported."),
    ] = "water",
    ph: Annotated[
        float,
        typer.Option("--ph", help="Candidate-generation pH recorded in metadata."),
    ] = 7.0,
) -> None:
    """Run CREST/xTB conformer search and ensemble reduction for seed conformers."""
    config = RunConfig(
        input_path=seeds_sdf,
        output_dir=out,
        chemistry={"ph": ph, "solvent": solvent},
    )
    seeds = read_seed_sdf(seeds_sdf)
    all_records = []
    try:
        for seed_record in seeds:
            all_records.extend(run_crest_for_seed(seed_record, config))
    except (CrestUnavailableError, CrestExecutionError) as exc:
        raise typer.BadParameter(str(exc)) from exc

    report = {
        "input_path": str(seeds_sdf),
        "seed_count": len(seeds),
        "crest_record_count": len(all_records),
        "failure_count": len([record for record in all_records if record.crest_index == 0]),
        "solvent": config.chemistry.solvent,
        "solvent_model": config.chemistry.solvent_model,
        "ph": config.chemistry.ph,
        "crest_ids": [record.id for record in all_records],
    }
    write_json(out / "crest" / "crest_report.json", report)
    console.print(f"CREST records: {len(all_records)}")
    console.print(f"Wrote outputs to [bold]{out / 'crest'}[/bold]")


@app.command(name="thermo")
def thermo_command(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
) -> None:
    """Run or inspect xTB thermochemistry records for a run directory."""
    config = RunConfig(input_path=run_dir, output_dir=run_dir)
    rankable_records = load_rankable_records(run_dir)
    if not rankable_records:
        raise typer.BadParameter(
            "No CREST conformer or xTB thermo records found. Run `dsvr run-crest` "
            "or the full `dsvr run` workflow before `dsvr thermo`."
        )
    if all(isinstance(record, ThermoRecord) for record in rankable_records):
        console.print(f"xTB thermo records already present: {len(rankable_records)}")
        console.print(f"Run directory: [bold]{run_dir}[/bold]")
        return

    thermo_records = []
    try:
        for record in rankable_records:
            if isinstance(record, CrestConformerRecord):
                thermo_records.append(run_xtb_thermo(record, config))
    except (XtbUnavailableError, XtbExecutionError) as exc:
        raise typer.BadParameter(
            f"xTB thermo failed: {exc}. Install xTB or run `dsvr doctor` for diagnostics."
        ) from exc

    write_json(
        run_dir / "xtb" / "thermo_report.json",
        {
            "run_dir": str(run_dir),
            "thermo_count": len(thermo_records),
            "thermo_ids": [record.id for record in thermo_records],
        },
    )
    console.print(f"xTB thermo records written: {len(thermo_records)}")
    console.print(f"Wrote outputs under [bold]{run_dir / 'xtb'}[/bold]")


@app.command(name="rank")
def rank_command(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
    population_scope: Annotated[
        str,
        typer.Option(
            "--population-scope",
            help="Population grouping: same_formula, same_charge, or all_approximate.",
        ),
    ] = "same_formula",
) -> None:
    """Compute ΔG values and approximate Boltzmann populations for a run directory."""
    config = RunConfig(
        input_path=run_dir,
        output_dir=run_dir,
        thermo={"population_scope": population_scope},
    )
    rankable_records = load_rankable_records(run_dir)
    if not rankable_records:
        raise typer.BadParameter(
            "No rankable CREST/xTB records found. Run `dsvr run-crest`, `dsvr thermo`, "
            "or the full `dsvr run` workflow first."
        )
    ranked = compute_delta_g_and_populations(rankable_records, config)
    output_dir = run_dir / "ranking"
    write_ranked_outputs(ranked, output_dir)
    write_final_ranked_outputs(run_dir, ranked, config)
    console.print(f"Ranked records: {len(ranked)}")
    console.print("Population scope: " + population_scope)
    console.print(
        "Population estimates are approximate across protomer/protonation states unless "
        "micro-pKa/proton chemical-potential corrections are available."
    )
    console.print(f"Wrote ranking outputs to [bold]{output_dir}[/bold]")
    console.print(f"Wrote final ranked outputs to [bold]{run_dir}[/bold]")


@app.command(name="run-qm")
def run_qm_command(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
    backend: Annotated[
        str,
        typer.Option("--backend", help="Optional QM backend: psi4, pyscf, or none."),
    ] = "none",
    top_n: Annotated[
        int,
        typer.Option("--top-n", help="Number of final ranked candidates to rescore."),
    ] = 5,
) -> None:
    """Optionally rescore the smallest final candidate set with Psi4 or PySCF."""
    config = RunConfig(
        input_path=run_dir,
        output_dir=run_dir,
        refinement={
            "qm_backend": backend,
            "psi4_enabled": backend == "psi4",
            "pyscf_enabled": backend == "pyscf",
            "max_candidates_for_refinement": top_n,
        },
    )
    records = load_ranked_for_qm(run_dir)
    if backend == "none":
        console.print("QM rescoring is optional and disabled with --backend none.")
        return
    try:
        if backend == "psi4":
            rescored = rescore_top_ranked_with_psi4(records, config)
        elif backend == "pyscf":
            rescored = rescore_top_ranked_with_pyscf(records, config)
        else:
            raise typer.BadParameter("backend must be psi4, pyscf, or none")
    except (Psi4UnavailableError, PySCFUnavailableError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    console.print(f"QM rescored records: {len(rescored)}")
    console.print(
        "QM rescoring is optional and uses electronic energies by default; prior "
        "rankings are preserved."
    )
    console.print(f"Wrote QM outputs to [bold]{run_dir / 'qm' / backend}[/bold]")


@app.command(name="summarize")
def summarize_command(
    run_dir: Annotated[Path, typer.Argument(exists=True, file_okay=False)],
    html: Annotated[
        bool,
        typer.Option("--html/--no-html", help="Also write a simple report.html file."),
    ] = False,
) -> None:
    """Regenerate a human-readable summary report for a run directory."""
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.exists():
        raise typer.BadParameter(
            f"{run_dir} does not contain manifest.json. Run `dsvr run` first."
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    resolved_config = run_dir / "resolved_config.yaml"
    config = (
        load_config(resolved_config)
        if resolved_config.exists()
        else RunConfig(input_path=run_dir, output_dir=run_dir)
    )
    ranked_records = _load_ranked_variants(run_dir)
    report_path = run_dir / "report.md"
    html_path = run_dir / "report.html" if html else None
    write_run_report(
        report_path,
        config=config,
        records=ranked_records,
        ranked_records=ranked_records,
        manifest=manifest,
        output_files=_known_output_files(run_dir),
        html_path=html_path,
    )
    console.print(f"Wrote summary report to [bold]{report_path}[/bold]")
    if html_path is not None:
        console.print(f"Wrote HTML report to [bold]{html_path}[/bold]")


def _warn_if_exhaustive(config: RunConfig, config_path: Path | None) -> None:
    if config.variant_filtering.mode == "exhaustive" or (
        config_path is not None and "exhaustive_debug" in config_path.name
    ):
        console.print(
            "[red]WARNING: exhaustive mode may generate extremely large numbers of "
            "variants and XYZ files. Use only for small molecules or debugging.[/red]"
        )


def _candidate_explosion_rows(config: RunConfig) -> list[tuple[str, str, int]]:
    protomers = config.protonation.max_protomers_per_molecule
    rdkit_tautomers = config.tautomer_filtering.max_rdkit_tautomers_before_auto3d
    selected_tautomers = config.tautomer_filtering.tauto_k
    stereoisomers = config.stereoisomer_filtering.max_stereoisomers_per_tautomer
    selected_stereoisomers = min(
        stereoisomers,
        config.stereoisomer_filtering.keep_top_n_diastereomers,
    )
    final_conformers = config.final_3d.k
    optional_validation = config.optional_validation.max_variants_per_molecule
    return [
        ("protomers", f"{protomers}", protomers),
        (
            "pre-Auto3D tautomer candidates",
            f"{protomers} protomers x {rdkit_tautomers} RDKit tautomers",
            protomers * rdkit_tautomers,
        ),
        (
            "selected tautomer states",
            f"{protomers} protomers x {selected_tautomers} tauto-k",
            protomers * selected_tautomers,
        ),
        (
            "stereoisomer candidates",
            f"{protomers * selected_tautomers} selected tautomers x {stereoisomers} stereo cap",
            protomers * selected_tautomers * stereoisomers,
        ),
        (
            "selected stereoisomer states",
            (
                f"{protomers * selected_tautomers} selected tautomers x "
                f"{selected_stereoisomers} stereo pruning cap"
            ),
            protomers * selected_tautomers * selected_stereoisomers,
        ),
        (
            "final Auto3D conformer jobs",
            (
                f"{protomers * selected_tautomers * selected_stereoisomers} variants x "
                f"{final_conformers} final k"
            ),
            protomers * selected_tautomers * selected_stereoisomers * final_conformers,
        ),
        (
            "optional CREST/xTB validation",
            f"enabled={config.optional_validation.crest_xtb_enabled}, cap={optional_validation}",
            optional_validation if config.optional_validation.crest_xtb_enabled else 0,
        ),
    ]


def _print_dry_run_summary(config: RunConfig, outdir: Path) -> None:
    console.print(f"Wrote dry-run plan to [bold]{outdir / 'dry_run_plan.json'}[/bold]")
    for step in planned_steps(config):
        tools = f" tools={','.join(step.external_tools)}" if step.external_tools else ""
        console.print(f"- {step.name}: {step.output_dir}{tools}")

    rows = _candidate_explosion_rows(config)
    console.print("Maximum candidate estimates per input molecule:")
    for label, formula, maximum in rows:
        console.print(f"- {label}: {formula} = {maximum}")

    table = Table(title="Maximum candidate estimates per input molecule")
    table.add_column("Stage cap")
    table.add_column("Formula")
    table.add_column("Maximum", justify="right")
    for label, formula, maximum in rows:
        table.add_row(label, formula, str(maximum))
    console.print(table)


def _execute_workflow_command(
    ctx: typer.Context,
    *,
    input_path: Path,
    config_path: Path | None,
    output_dir: Path | None,
    dry_run: bool,
    ph: float | None,
    solvent: str | None,
    max_protomers: int | None = None,
    tauto_k: int | None = None,
    tauto_window: float | None = None,
    max_stereoisomers: int | None = None,
    seeding_method: SeederMethod | None = None,
    censo_enabled: bool | None = None,
    crest_xtb_enabled: bool | None = None,
    agent_enabled: bool | None = None,
    overwrite: bool | None = None,
    resume: bool | None = None,
    force_ligprep_like: bool = False,
) -> None:
    globals_ = _global_options(ctx)
    effective_config_path = config_path or globals_.get("config_path")
    effective_dry_run = dry_run or bool(globals_.get("dry_run"))
    config = load_config(effective_config_path) if effective_config_path else RunConfig()
    config = merge_cli_overrides(
        config,
        input_path=input_path,
        output_dir=output_dir,
        workflow_mode="ligprep_like" if force_ligprep_like else None,
        ph=ph,
        solvent=solvent,
        max_protomers=max_protomers,
        tauto_k=tauto_k,
        tauto_window=tauto_window,
        max_stereoisomers=max_stereoisomers,
        seeding_method=seeding_method,
        censo_enabled=censo_enabled,
        crest_xtb_enabled=crest_xtb_enabled,
        agent_enabled=agent_enabled,
    )
    config = _apply_runtime_options(
        config,
        globals_=globals_,
        dry_run=effective_dry_run,
        overwrite=overwrite,
        resume=resume,
    )
    _warn_if_exhaustive(config, effective_config_path)
    result = run_workflow(config=config)
    if config.dry_run:
        _print_dry_run_summary(config, result.outdir)
    else:
        console.print(f"Wrote workflow outputs to [bold]{result.outdir}[/bold]")
    console.print(f"Molecules: {result.molecule_count}")


@app.command(name="prepare-ligands")
def prepare_ligands(
    ctx: typer.Context,
    input_path: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    config_path: Annotated[
        Path | None,
        typer.Option(
            "--config",
            "-c",
            exists=True,
            dir_okay=False,
            help="YAML ligand-prep config.",
        ),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option("--out", "--outdir", "-o", help="Run output directory."),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            help="Show planned ligand-prep stages and maximum candidate counts.",
        ),
    ] = False,
    ph: Annotated[
        float | None,
        typer.Option("--ph", help="Override ligand-preparation pH."),
    ] = None,
    solvent: Annotated[
        str | None,
        typer.Option("--solvent", help="Override solvent name passed to configured tools."),
    ] = None,
    max_protomers: Annotated[
        int | None,
        typer.Option("--max-protomers", help="Cap plausible protomers per input molecule."),
    ] = None,
    tauto_k: Annotated[
        int | None,
        typer.Option("--tauto-k", help="Keep top K Auto3D-ranked tautomers per protomer."),
    ] = None,
    tauto_window: Annotated[
        float | None,
        typer.Option(
            "--tauto-window",
            help="Keep tautomers within this Auto3D energy window in kcal/mol.",
        ),
    ] = None,
    max_stereoisomers: Annotated[
        int | None,
        typer.Option("--max-stereoisomers", help="Cap stereoisomers per selected tautomer."),
    ] = None,
    enable_crest_validation: Annotated[
        bool,
        typer.Option(
            "--enable-crest-validation",
            help="Run optional CREST/xTB validation on selected final variants.",
        ),
    ] = False,
    agent: Annotated[
        bool,
        typer.Option(
            "--agent",
            help="Enable experimental local agent diagnostics only on failures.",
        ),
    ] = False,
    overwrite: Annotated[
        bool | None,
        typer.Option("--overwrite/--no-overwrite", help="Override workflow overwrite policy."),
    ] = None,
    resume: Annotated[
        bool | None,
        typer.Option("--resume/--no-resume", help="Override workflow resume policy."),
    ] = None,
) -> None:
    """Prepare ligand structural variants with the default LigPrep-like workflow."""
    _execute_workflow_command(
        ctx,
        input_path=input_path,
        config_path=config_path,
        output_dir=output_dir,
        dry_run=dry_run,
        ph=ph,
        solvent=solvent,
        max_protomers=max_protomers,
        tauto_k=tauto_k,
        tauto_window=tauto_window,
        max_stereoisomers=max_stereoisomers,
        crest_xtb_enabled=True if enable_crest_validation else None,
        agent_enabled=True if agent else None,
        overwrite=overwrite,
        resume=resume,
        force_ligprep_like=True,
    )


@app.command()
def run(
    ctx: typer.Context,
    input_path: Annotated[Path, typer.Argument(exists=True, dir_okay=False)],
    config_path: Annotated[
        Path | None,
        typer.Option("--config", "-c", exists=True, dir_okay=False, help="YAML workflow config."),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option("--out", "--outdir", "-o", help="Run output directory."),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Show planned workflow steps without external execution."),
    ] = False,
    ph: Annotated[
        float | None,
        typer.Option("--ph", help="Override candidate-generation pH."),
    ] = None,
    solvent: Annotated[
        str | None,
        typer.Option("--solvent", help="Override solvent name passed to configured tools."),
    ] = None,
    max_protomers: Annotated[
        int | None,
        typer.Option("--max-protomers", help="Cap plausible protomers per input molecule."),
    ] = None,
    tauto_k: Annotated[
        int | None,
        typer.Option("--tauto-k", help="Keep top K Auto3D-ranked tautomers per protomer."),
    ] = None,
    tauto_window: Annotated[
        float | None,
        typer.Option(
            "--tauto-window",
            help="Keep tautomers within this Auto3D energy window in kcal/mol.",
        ),
    ] = None,
    max_stereoisomers: Annotated[
        int | None,
        typer.Option("--max-stereoisomers", help="Cap stereoisomers per selected tautomer."),
    ] = None,
    seeding_method: Annotated[
        SeederMethod | None,
        typer.Option("--seeding-method", help="Override seeding method."),
    ] = None,
    censo_enabled: Annotated[
        bool | None,
        typer.Option("--enable-censo/--disable-censo", help="Override optional CENSO refinement."),
    ] = None,
    crest_xtb_enabled: Annotated[
        bool | None,
        typer.Option(
            "--enable-crest-validation/--disable-crest-validation",
            help="Run optional CREST/xTB validation on selected final variants.",
        ),
    ] = None,
    agent: Annotated[
        bool,
        typer.Option(
            "--agent",
            help="Enable experimental local agent diagnostics only on failures.",
        ),
    ] = False,
    overwrite: Annotated[
        bool | None,
        typer.Option("--overwrite/--no-overwrite", help="Override workflow overwrite policy."),
    ] = None,
    resume: Annotated[
        bool | None,
        typer.Option("--resume/--no-resume", help="Override workflow resume policy."),
    ] = None,
) -> None:
    """Run the full DSVR orchestration workflow; prefer prepare-ligands for ligand prep."""
    _execute_workflow_command(
        ctx,
        input_path=input_path,
        config_path=config_path,
        output_dir=output_dir,
        dry_run=dry_run,
        ph=ph,
        solvent=solvent,
        max_protomers=max_protomers,
        tauto_k=tauto_k,
        tauto_window=tauto_window,
        max_stereoisomers=max_stereoisomers,
        seeding_method=seeding_method,
        censo_enabled=censo_enabled,
        crest_xtb_enabled=crest_xtb_enabled,
        agent_enabled=True if agent else None,
        overwrite=overwrite,
        resume=resume,
    )


def _load_ranked_variants(run_dir: Path) -> list[RankedVariantRecord]:
    for path in (run_dir / "ranked_variants.json", run_dir / "ranking" / "ranked_variants.json"):
        if path.exists():
            payload = json.loads(path.read_text(encoding="utf-8"))
            return [RankedVariantRecord.model_validate(item) for item in payload]
    return []


def _known_output_files(run_dir: Path) -> list[Path]:
    names = [
        "manifest.json",
        "resolved_config.yaml",
        "logs",
        "invalid_inputs.csv",
        "inputs.csv",
        "protomers.csv",
        "tautomers.csv",
        "stereoisomers.csv",
        "seeds.csv",
        "crest_conformers.csv",
        "thermo.csv",
        "ranked_variants.csv",
        "ranked_variants.json",
        "ranked_variants.sdf",
        "report.md",
        "report.html",
    ]
    return [run_dir / name for name in names if (run_dir / name).exists()]


if __name__ == "__main__":
    app()
