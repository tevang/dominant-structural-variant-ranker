from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from dsvr.reporting.progress import read_jsonl_tail

ACTIVE_STATUSES = {"started", "running", "progress"}
DONE_STATUSES = {"completed", "skipped"}


def run_status(run_dir: Path, *, log_tail_lines: int = 20) -> dict[str, Any]:
    progress = _read_json(run_dir / "progress.json")
    last_event = progress.get("last_event", {}) if isinstance(progress, dict) else {}
    stage_rows = _stage_rows_from_csv(run_dir / "stage_summary.csv")
    stage_counts = progress.get("stage_counts", {}) if isinstance(progress, dict) else {}
    active_stage = last_event.get("stage") if last_event.get("status") in ACTIVE_STATUSES else None
    last_completed_stage = _last_completed_stage(stage_rows)
    warnings = read_jsonl_tail(run_dir / "warnings.jsonl", limit=5)
    failures = read_jsonl_tail(run_dir / "failures.jsonl", limit=5)
    recovery_failures = read_jsonl_tail(run_dir / "checkpoints" / "failures.jsonl", limit=5)
    molecule_states = _molecule_state_rows(run_dir / "checkpoints" / "molecules")

    return {
        "run_dir": str(run_dir),
        "last_stage": last_event.get("stage"),
        "last_status": last_event.get("status"),
        "last_completed_stage": last_completed_stage,
        "active_stage": active_stage,
        "current_molecule": last_event.get("molecule_name"),
        "current_molecule_index": last_event.get("molecule_index"),
        "current_molecule_total": last_event.get("molecule_total"),
        "last_completed_molecule": last_event.get("molecule_name"),
        "active_command": last_event.get("active_command"),
        "disk_usage_mb": _directory_size_mb(run_dir),
        "xyz_file_count": _xyz_file_count(run_dir),
        "stage_counts": stage_counts or {row["stage"]: 1 for row in stage_rows if row.get("stage")},
        "counts_by_stage": stage_rows,
        "latest_warnings": warnings,
        "latest_failures": failures,
        "latest_recovery_failures": recovery_failures,
        "molecule_state_count": len(molecule_states),
        "molecule_states": molecule_states[-10:],
        "warning_count": (
            progress.get("warning_count", len(warnings))
            if isinstance(progress, dict)
            else len(warnings)
        ),
        "failure_count": (
            progress.get("failure_count", len(failures))
            if isinstance(progress, dict)
            else len(failures)
        ),
        "resume_possible": _resume_possible(run_dir, last_event),
        "latest_log_tail": _latest_log_tail(run_dir, log_tail_lines),
    }


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _stage_rows_from_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle) if row.get("stage")]


def _last_completed_stage(rows: list[dict[str, str]]) -> str | None:
    completed = [row.get("stage") for row in rows if row.get("status") in DONE_STATUSES]
    return completed[-1] if completed else None


def _resume_possible(run_dir: Path, last_event: dict[str, Any]) -> bool:
    if not run_dir.exists():
        return False
    if last_event.get("stage") == "Report writing" and last_event.get("status") == "completed":
        return False
    return (run_dir / "resolved_config.yaml").exists() or any(run_dir.glob("**/done.json"))


def _molecule_state_rows(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
    for item in sorted(path.glob("*.json")):
        try:
            rows.append(json.loads(item.read_text(encoding="utf-8")))
        except json.JSONDecodeError:
            continue
    return rows


def _latest_log_tail(run_dir: Path, lines: int) -> str:
    logs = sorted(run_dir.glob("**/*.log"), key=lambda path: path.stat().st_mtime)
    if not logs:
        return ""
    latest = logs[-1]
    text = latest.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(text[-lines:])


def _directory_size_mb(path: Path) -> float:
    if not path.exists():
        return 0.0
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file():
                total += item.stat().st_size
        except OSError:
            continue
    return total / (1024**2)


def _xyz_file_count(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for item in path.rglob("*.xyz") if item.is_file())
