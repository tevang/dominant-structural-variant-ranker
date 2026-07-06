from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from dsvr.config import RunConfig


@dataclass(frozen=True)
class BugPackage:
    path: Path
    failure: dict[str, Any]
    command: dict[str, Any]
    prompt_context: str


def build_bug_package(run_dir: Path, config: RunConfig, *, latest: bool = True) -> BugPackage:
    if not latest:
        raise ValueError("Only --latest bug package selection is currently supported")
    package_dir = run_dir / "bug_package"
    package_dir.mkdir(parents=True, exist_ok=True)

    failure = _latest_jsonl_record(run_dir / "checkpoints" / "failures.jsonl")
    if not failure:
        failure = _latest_jsonl_record(run_dir / "failures.jsonl")
    if not failure:
        failure = {
            "failure_kind": "UNKNOWN",
            "stage": "",
            "item_id": None,
            "item_name": None,
            "message": "No recorded failure found in run directory.",
        }
    _write_json(package_dir / "failure.json", failure)

    command = _latest_command_record(run_dir)
    _write_json(package_dir / "command.json", command)
    (package_dir / "last_200_lines.log").write_text(
        _last_200_log_lines(run_dir, command),
        encoding="utf-8",
    )
    (package_dir / "config_fragment.yaml").write_text(
        yaml.safe_dump(_config_fragment(config), sort_keys=False),
        encoding="utf-8",
    )
    (package_dir / "molecule.smi_or_sdf").write_text(
        _molecule_fragment(failure, config),
        encoding="utf-8",
    )
    _copy_stage_summary(run_dir, package_dir)

    context = _compact_context(package_dir)
    return BugPackage(path=package_dir, failure=failure, command=command, prompt_context=context)


def _latest_jsonl_record(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    latest: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(payload, dict):
            latest = payload
    return latest


def _latest_command_record(run_dir: Path) -> dict[str, Any]:
    candidates = sorted(
        run_dir.glob("**/command.json"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    for path in candidates:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(payload, dict):
            payload.setdefault("command_json", str(path))
            return payload
    return {}


def _last_200_log_lines(run_dir: Path, command: dict[str, Any]) -> str:
    log_paths = [
        Path(str(command[key]))
        for key in ("stderr_log", "stdout_log")
        if command.get(key)
    ]
    if command.get("command_json"):
        log_paths.append(Path(str(command["command_json"])).with_name("combined.log"))
    if not log_paths:
        log_paths = sorted(
            run_dir.glob("**/*.log"),
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )[:1]
    lines: list[str] = []
    for path in log_paths:
        if not path.exists() or not path.is_file():
            continue
        text = path.read_text(encoding="utf-8", errors="replace").splitlines()
        lines.extend([f"== {path.name} =="] + text[-200:])
    return "\n".join(lines[-200:]) + ("\n" if lines else "")


def _config_fragment(config: RunConfig) -> dict[str, Any]:
    data = config.model_dump(mode="json")
    keys = [
        "agent",
        "error_handling",
        "enumeration",
        "tautomer_filtering",
        "stereoisomer_filtering",
        "seeding",
        "final_3d",
        "logging",
    ]
    return {key: data[key] for key in keys if key in data}


def _molecule_fragment(failure: dict[str, Any], config: RunConfig) -> str:
    lines = [
        f"item_id: {failure.get('item_id') or ''}",
        f"item_name: {failure.get('item_name') or ''}",
        f"input_path: {config.input_path}",
    ]
    try:
        input_lines = config.input_path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()
    except OSError:
        input_lines = []
    if input_lines:
        lines.append("input_head:")
        lines.extend(input_lines[:20])
    return "\n".join(lines) + "\n"


def _copy_stage_summary(run_dir: Path, package_dir: Path) -> None:
    source = run_dir / "stage_summary.csv"
    target = package_dir / "stage_summary.csv"
    if source.exists():
        shutil.copyfile(source, target)
    else:
        target.write_text("", encoding="utf-8")


def _compact_context(package_dir: Path) -> str:
    parts = []
    for filename in (
        "failure.json",
        "command.json",
        "last_200_lines.log",
        "config_fragment.yaml",
        "molecule.smi_or_sdf",
        "stage_summary.csv",
    ):
        path = package_dir / filename
        text = path.read_text(encoding="utf-8", errors="replace")
        parts.append(f"## {filename}\n{text[:4000]}")
    return "\n\n".join(parts)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
