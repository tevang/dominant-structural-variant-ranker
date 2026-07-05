from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from dsvr.config import RunConfig
from dsvr.utils.hashing import sha256_text


@dataclass(frozen=True)
class WorkflowStep:
    name: str
    description: str
    output_dir: Path
    external_tools: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class StepState:
    step: WorkflowStep
    status: str
    skipped: bool
    output_dir: Path
    input_hash: str
    config_hash: str
    details: dict[str, Any] = field(default_factory=dict)


WORKFLOW_STEP_NAMES = [
    "input",
    "standardize",
    "protonation",
    "tautomers",
    "stereochemistry",
    "seeding",
    "crest",
    "xtb_thermo",
    "ranking",
    "censo",
    "qm",
    "reports",
]


def planned_steps(config: RunConfig) -> list[WorkflowStep]:
    root = config.output_dir
    return [
        WorkflowStep("input", "Read and validate input molecules.", root / "input"),
        WorkflowStep("standardize", "Standardize RDKit molecules when enabled.", root / "input"),
        WorkflowStep(
            "protonation",
            "Generate pH/protomer candidates with molscrub.",
            root / "enumeration" / "protomers",
            ["molscrub"],
        ),
        WorkflowStep(
            "tautomers",
            (
                "Auto3D tautomer enumeration."
                if config.protocol == "auto3d_entropy"
                else "Enumerate RDKit tautomer candidates."
            ),
            root / "enumeration" / "tautomers",
            ["Auto3D"] if config.protocol == "auto3d_entropy" else [],
        ),
        WorkflowStep(
            "stereochemistry",
            (
                "Auto3D stereoisomer enumeration."
                if config.protocol == "auto3d_entropy"
                else "Enumerate explicit RDKit stereoisomer candidates."
            ),
            root / "enumeration" / "stereoisomers",
            ["Auto3D"] if config.protocol == "auto3d_entropy" else [],
        ),
        WorkflowStep(
            "seeding",
            (
                "Run Auto3D representative conformer generation."
                if config.protocol == "auto3d_entropy"
                else f"Generate 3D seeds with {config.seeding.method}."
            ),
            root / "seeding",
            ["Auto3D"]
            if config.protocol == "auto3d_entropy" or config.seeding.method in {"auto3d", "both"}
            else [],
        ),
        WorkflowStep(
            "crest",
            (
                "Skipped by Auto3D entropy protocol."
                if config.protocol == "auto3d_entropy"
                else "Run CREST/xTB conformer search and ensemble reduction."
            ),
            root / "crest",
            [] if config.protocol == "auto3d_entropy" else ["crest", "xtb"],
        ),
        WorkflowStep(
            "xtb_thermo",
            (
                "Skipped by Auto3D representative protocol."
                if config.protocol == "auto3d_entropy"
                else "Run xTB Hessian/thermo and collect free-energy estimates."
            ),
            root / "xtb",
            [] if config.protocol == "auto3d_entropy" else ["xtb"],
        ),
        WorkflowStep("ranking", "Rank variants and approximate populations.", root / "ranking"),
        WorkflowStep(
            "censo",
            "Optionally refine top candidates with CENSO.",
            root / "censo",
            ["censo"],
        ),
        WorkflowStep(
            "qm",
            "Optionally rescore final candidates with Psi4/PySCF.",
            root / "qm",
            [config.refinement.qm_backend]
            if config.refinement.qm_backend in {"psi4", "pyscf"}
            else [],
        ),
        WorkflowStep("reports", "Write manifest and summary reports.", root),
    ]


def config_hash(config: RunConfig) -> str:
    return sha256_text(json.dumps(config.model_dump(mode="json"), sort_keys=True, default=str))


def file_hash(path: Path) -> str:
    return sha256_text(path.read_text(encoding="utf-8", errors="replace"))


def records_hash(records: list[object]) -> str:
    payload = [
        getattr(record, "id", None)
        or getattr(record, "input_id", None)
        or getattr(record, "molname", None)
        or repr(record)
        for record in records
    ]
    return sha256_text(json.dumps(payload, sort_keys=True, default=str))


def should_skip_step(step: WorkflowStep, input_hash: str, config: RunConfig) -> bool:
    if config.overwrite or not config.resume:
        return False
    done = read_done(step.output_dir)
    return (
        done.get("step") == step.name
        and done.get("input_hash") == input_hash
        and done.get("config_hash") == config_hash(config)
        and _step_artifacts_complete(step, done, config)
    )


def mark_done(
    step: WorkflowStep,
    input_hash: str,
    config: RunConfig,
    *,
    details: dict[str, Any] | None = None,
) -> StepState:
    step.output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "step": step.name,
        "status": "done",
        "input_hash": input_hash,
        "config_hash": config_hash(config),
        "completed_at": datetime.now(timezone.utc).isoformat(),  # noqa: UP017
        "details": details or {},
    }
    (step.output_dir / "done.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    return StepState(
        step=step,
        status="done",
        skipped=False,
        output_dir=step.output_dir,
        input_hash=input_hash,
        config_hash=payload["config_hash"],
        details=details or {},
    )


def skipped_state(step: WorkflowStep, input_hash: str, config: RunConfig) -> StepState:
    return StepState(
        step=step,
        status="skipped",
        skipped=True,
        output_dir=step.output_dir,
        input_hash=input_hash,
        config_hash=config_hash(config),
        details=read_done(step.output_dir).get("details", {}),
    )


def read_done(output_dir: Path) -> dict[str, Any]:
    path = output_dir / "done.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _step_artifacts_complete(
    step: WorkflowStep,
    done: dict[str, Any],
    config: RunConfig,
) -> bool:
    details = done.get("details", {})
    count = details.get("count")
    if step.name == "standardize":
        return (step.output_dir / "standardized_inputs.csv").exists()
    if step.name == "protonation":
        return count != 0 and any(step.output_dir.glob("*_protomers.sdf"))
    if step.name == "tautomers":
        return count != 0 and any(step.output_dir.glob("*_tautomers.sdf"))
    if step.name == "stereochemistry":
        return count != 0 and any(step.output_dir.glob("*_stereoisomers.sdf"))
    if step.name == "seeding":
        return count != 0 and any(step.output_dir.glob("**/*_seeds.sdf"))
    if step.name == "crest":
        if not config.crest.enabled:
            return False
        return count != 0 and _has_successful_crest_records(step.output_dir)
    if step.name == "xtb_thermo":
        input_count = details.get("input_count", 0)
        if not config.thermo.enabled or input_count == 0:
            return True
        return count != 0 and any(step.output_dir.glob("**/xtb_thermo.json"))
    if step.name == "ranking":
        return count != 0 and (step.output_dir / "ranked_variants.json").exists()
    if step.name == "reports":
        return (step.output_dir / "manifest.json").exists() and (
            step.output_dir / "summary.md"
        ).exists()
    return True


def _has_successful_crest_records(output_dir: Path) -> bool:
    for path in output_dir.glob("**/crest_provenance.jsonl"):
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            if record.get("crest_index", 0) and record.get("energy_kcal_mol") is not None:
                return True
    return False


def write_dry_run_plan(config: RunConfig) -> Path:
    path = config.output_dir / "dry_run_plan.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = [
        {
            "name": step.name,
            "description": step.description,
            "output_dir": str(step.output_dir),
            "external_tools": step.external_tools,
        }
        for step in planned_steps(config)
    ]
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path
