from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


def run_status(run_dir: Path, *, log_tail_lines: int = 20) -> dict[str, Any]:
    progress = _read_json(run_dir / "progress.json")
    last_event = progress.get("last_event", {}) if isinstance(progress, dict) else {}
    stage_counts = progress.get("stage_counts", {}) if isinstance(progress, dict) else {}
    return {
        "run_dir": str(run_dir),
        "last_stage": last_event.get("stage"),
        "last_status": last_event.get("status"),
        "last_completed_molecule": last_event.get("molecule_name"),
        "active_command": last_event.get("active_command"),
        "disk_usage_mb": _directory_size_mb(run_dir),
        "xyz_file_count": _xyz_file_count(run_dir),
        "stage_counts": stage_counts or _stage_counts_from_csv(run_dir / "stage_summary.csv"),
        "latest_log_tail": _latest_log_tail(run_dir, log_tail_lines),
    }


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _stage_counts_from_csv(path: Path) -> dict[str, int]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        return {row["stage"]: 1 for row in csv.DictReader(handle) if row.get("stage")}


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
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file()) / (1024**2)


def _xyz_file_count(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for item in path.rglob("*.xyz") if item.is_file())
