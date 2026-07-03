from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path
from typing import Any

from dsvr.models import ToolStatus
from dsvr.runners.subprocess_utils import (
    executable_version,
    meets_minimum_version,
    python_import_check,
    which_executable,
)

PYTHON_MINIMUM = "3.11"
MODULE_CHECKS = {
    "rdkit": {"required": True, "minimum": None},
    "molscrub": {"required": True, "minimum": None},
    "Auto3D": {"required": False, "minimum": None},
    "psi4": {"required": False, "minimum": None},
    "pyscf": {"required": False, "minimum": None},
}
EXECUTABLE_CHECKS = {
    "scrub.py": {"required": False, "version_args": ["-h"]},
    "molscrub": {"required": False, "version_args": ["-h"]},
    "auto3d": {"required": False, "version_args": ["--version"]},
    "auto3D": {"required": False, "version_args": ["--version"]},
    "xtb": {"required": True, "version_args": ["--version"]},
    "crest": {"required": True, "version_args": ["--version"]},
    "censo": {"required": False, "version_args": ["--version"]},
    "psi4": {"required": False, "version_args": ["--version"]},
}


def check_tools(output_dir: Path | None = None) -> list[ToolStatus]:
    statuses = [
        _python_status(),
        *_module_statuses(),
        *_executable_statuses(),
        *_compound_tool_statuses(),
        _writable_output_status(output_dir or Path("runs/dsvr")),
        _cpu_status(),
        _disk_status(output_dir or Path.cwd()),
    ]
    return statuses


def doctor_payload(output_dir: Path | None = None) -> dict[str, Any]:
    statuses = check_tools(output_dir=output_dir)
    required_missing = [
        status.name for status in statuses if status.required and not status.available
    ]
    return {
        "ok": not required_missing,
        "strict_failure_count": len(required_missing),
        "required_missing": required_missing,
        "checks": [status.model_dump(mode="json") for status in statuses],
    }


def _python_status() -> ToolStatus:
    version = ".".join(str(part) for part in sys.version_info[:3])
    minimum_ok = meets_minimum_version(version, PYTHON_MINIMUM)
    return ToolStatus(
        name="python",
        kind="runtime",
        required=True,
        available=bool(minimum_ok),
        detail=sys.executable,
        version=version,
        minimum_version=PYTHON_MINIMUM,
        meets_minimum_version=minimum_ok,
    )


def _module_statuses() -> list[ToolStatus]:
    statuses: list[ToolStatus] = []
    for module_name, config in MODULE_CHECKS.items():
        available, version = python_import_check(module_name)
        minimum = config["minimum"]
        minimum_ok = meets_minimum_version(version, minimum)
        detail = "importable" if available else _module_install_hint(module_name)
        statuses.append(
            ToolStatus(
                name=module_name,
                kind="python-module",
                required=bool(config["required"]),
                available=available and minimum_ok is not False,
                detail=detail,
                version=version,
                minimum_version=minimum,
                meets_minimum_version=minimum_ok,
            )
        )
    return statuses


def _executable_statuses() -> list[ToolStatus]:
    statuses: list[ToolStatus] = []
    for executable, config in EXECUTABLE_CHECKS.items():
        path = which_executable(executable)
        version = None
        if path is not None:
            version = executable_version(executable, args=config["version_args"])
        statuses.append(
            ToolStatus(
                name=executable,
                kind="executable",
                required=bool(config["required"]),
                available=path is not None,
                detail=path or "not on PATH",
                version=version,
            )
        )
    return statuses


def _compound_tool_statuses() -> list[ToolStatus]:
    molscrub_module, molscrub_version = python_import_check("molscrub")
    molscrub_cli = which_executable("scrub.py") or which_executable("molscrub")
    auto3d_module, auto3d_version = python_import_check("Auto3D")
    auto3d_cli = (
        which_executable("auto3d")
        or which_executable("auto3D")
        or which_executable("Auto3D")
    )
    return [
        ToolStatus(
            name="molscrub-api-or-cli",
            kind="compound",
            required=True,
            available=molscrub_module or molscrub_cli is not None,
            detail=(
                "Python import or CLI available"
                if molscrub_module or molscrub_cli is not None
                else "install molscrub Python package or provide scrub.py/molscrub on PATH"
            ),
            version=molscrub_version,
        ),
        ToolStatus(
            name="Auto3D-api-or-cli",
            kind="compound",
            required=False,
            available=auto3d_module or auto3d_cli is not None,
            detail=(
                "Python import or CLI available"
                if auto3d_module or auto3d_cli is not None
                else "optional; install Auto3D Python package or provide auto3d on PATH"
            ),
            version=auto3d_version,
        ),
    ]


def _writable_output_status(output_dir: Path) -> ToolStatus:
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        probe = output_dir / ".dsvr_write_test"
        probe.write_text("ok\n", encoding="utf-8")
        probe.unlink()
        return ToolStatus(
            name="output-directory",
            kind="filesystem",
            required=True,
            available=True,
            detail=str(output_dir),
        )
    except OSError as exc:
        return ToolStatus(
            name="output-directory",
            kind="filesystem",
            required=True,
            available=False,
            detail=f"{output_dir}: {type(exc).__name__}: {exc}",
        )


def _cpu_status() -> ToolStatus:
    count = os.cpu_count() or 0
    return ToolStatus(
        name="cpu-count",
        kind="system",
        required=False,
        available=count > 0,
        detail=str(count),
    )


def _disk_status(path: Path) -> ToolStatus:
    target = path if path.exists() else path.parent
    usage = shutil.disk_usage(target)
    free_gb = usage.free / 1024**3
    return ToolStatus(
        name="disk-space",
        kind="system",
        required=False,
        available=usage.free > 0,
        detail=f"{free_gb:.2f} GiB free at {target}",
    )


def _module_install_hint(module_name: str) -> str:
    if module_name == "molscrub":
        return (
            "not importable; install with "
            "pip install git+https://github.com/forlilab/molscrub.git"
        )
    if module_name == "Auto3D":
        return "not importable; install with pip install Auto3D"
    return "not importable"
