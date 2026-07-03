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
            "Enumerate RDKit tautomer candidates.",
            root / "enumeration" / "tautomers",
        ),
        WorkflowStep(
            "stereochemistry",
            "Enumerate explicit RDKit stereoisomer candidates.",
            root / "enumeration" / "stereoisomers",
        ),
        WorkflowStep(
            "seeding",
            f"Generate 3D seeds with {config.seeding.method}.",
            root / "seeding",
            ["Auto3D"] if config.seeding.method in {"auto3d", "both"} else [],
        ),
        WorkflowStep(
            "crest",
            "Run CREST/xTB conformer search and ensemble reduction.",
            root / "crest",
            ["crest", "xtb"],
        ),
        WorkflowStep(
            "xtb_thermo",
            "Run xTB Hessian/thermo and collect free-energy estimates.",
            root / "xtb",
            ["xtb"],
        ),
        WorkflowStep("ranking", "Compute ΔG and approximate populations.", root / "ranking"),
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
