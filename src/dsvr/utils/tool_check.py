from __future__ import annotations

import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover
    from dsvr.config import ProtonationConfig

from dsvr.models import ToolStatus
from dsvr.runners.subprocess_utils import (
    executable_version,
    meets_minimum_version,
    python_import_check,
    which_executable,
)

PYTHON_MINIMUM = "3.11"

# Single-interface tools: (name, required).
_MODULE_CHECKS = {
    "rdkit": {"required": True},
    "pyscf": {"required": False},
}
_EXECUTABLE_CHECKS = {
    "xtb": {"required": True, "version_args": ["--version"]},
    "crest": {"required": True, "version_args": ["--version"]},
    "censo": {"required": False, "version_args": ["--version"]},
}


@dataclass(frozen=True)
class _InterfaceCheck:
    """One way to access a tool: a Python module import or a CLI executable."""

    kind: str  # "python-module" or "executable"
    label: str  # human-readable interface label, e.g. "python module" or "CLI"
    probe_names: tuple[str, ...]  # module name, or executable candidates in PATH order
    version_args: tuple[str, ...] = ("--version",)


@dataclass(frozen=True)
class _ToolGroup:
    """A logical external tool that is usable if any of its interfaces works."""

    name: str
    required: bool
    interfaces: tuple[_InterfaceCheck, ...]
    install_hint: str


_TOOL_GROUPS = (
    _ToolGroup(
        name="molscrub",
        required=False,
        interfaces=(
            _InterfaceCheck(
                kind="python-module",
                label="python module",
                probe_names=("molscrub",),
            ),
            _InterfaceCheck(
                kind="executable",
                label="CLI",
                probe_names=("scrub.py", "molscrub"),
                version_args=("-h",),
            ),
        ),
        install_hint=(
            "optional (legacy protonation tool; only needed when "
            "protonation.tool=molscrub); install with "
            "pip install git+https://github.com/forlilab/molscrub.git "
            "or provide scrub.py/molscrub on PATH"
        ),
    ),
    _ToolGroup(
        name="Auto3D",
        required=False,
        interfaces=(
            _InterfaceCheck(
                kind="python-module",
                label="python module",
                probe_names=("Auto3D",),
            ),
            _InterfaceCheck(
                kind="executable",
                label="CLI",
                probe_names=("auto3d", "auto3D", "Auto3D"),
            ),
        ),
        install_hint="optional; install with pip install Auto3D or provide auto3d on PATH",
    ),
    _ToolGroup(
        name="psi4",
        required=False,
        interfaces=(
            _InterfaceCheck(
                kind="python-module",
                label="python module",
                probe_names=("psi4",),
            ),
            _InterfaceCheck(
                kind="executable",
                label="CLI",
                probe_names=("psi4",),
            ),
        ),
        install_hint=(
            "optional; install the psi4 Python package or provide the psi4 "
            "executable on PATH"
        ),
    ),
)


def check_tools(
    output_dir: Path | None = None,
    protonation: ProtonationConfig | None = None,
) -> list[ToolStatus]:
    """Tool availability rows; `protonation` reflects a user config when given.

    Without a config, the default protonation selection (unipka, enabled) is
    assumed, which keeps Uni-Pka a required check.
    """

    statuses = [
        _python_status(),
        *_module_statuses(),
        _unipka_status(protonation),
        *_tool_group_statuses(),
        *_executable_statuses(),
        _writable_output_status(output_dir or Path("runs/dsvr")),
        _cpu_status(),
        _disk_status(output_dir or Path.cwd()),
    ]
    return statuses


def _unipka_status(protonation: ProtonationConfig | None = None) -> ToolStatus:
    """Uni-Pka (EasyDock container) availability: runtime + configured image.

    Required only when the (optionally user-selected) protonation stage uses
    Uni-Pka; with molscrub selected or protonation disabled, Uni-Pka becomes an
    informational row so `dsvr doctor --strict` reflects the actual workflow.
    """

    from dsvr.config import ProtonationConfig as _ProtonationConfig
    from dsvr.runners.unipka_runner import inspect_unipka

    config = protonation or _ProtonationConfig()
    selected = config.enabled and config.tool == "unipka"
    if not selected:
        reason = (
            "protonation disabled"
            if not config.enabled
            else f"protonation.tool={config.tool}; Uni-Pka not used"
        )
        return ToolStatus(
            name="unipka",
            kind="tool",
            required=False,
            available=True,
            detail=reason,
        )
    probe = inspect_unipka(config.unipka)
    runtime_ok = probe["runtime"] is not None
    image = probe["image"]
    available = runtime_ok and image is not None
    if available:
        detail = f"container runtime={probe['runtime']}; image={image}"
        script_override = probe.get("script_override")
        if script_override:
            detail += f"; script={script_override}"
    elif not runtime_ok:
        detail = str(probe["runtime_error"])
    else:
        detail = (
            f"Uni-Pka image not found: {probe['container']}; download unipka.sif "
            "from Zenodo or build from the EasyDock recipe (docs/external_tools.md)"
        )
    return ToolStatus(
        name="unipka",
        kind="tool",
        required=True,
        available=available,
        detail=detail,
    )


def doctor_payload(
    output_dir: Path | None = None,
    protonation: ProtonationConfig | None = None,
) -> dict[str, Any]:
    statuses = check_tools(output_dir=output_dir, protonation=protonation)
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
    for module_name, config in _MODULE_CHECKS.items():
        available, version = python_import_check(module_name)
        statuses.append(
            ToolStatus(
                name=module_name,
                kind="python-module",
                required=bool(config["required"]),
                available=available,
                detail="importable" if available else "not importable",
                version=version,
            )
        )
    return statuses


def _tool_group_statuses() -> list[ToolStatus]:
    """One summary row per tool followed by one row per interface.

    The summary row reports whether the tool is usable at all (any interface
    available); interface rows are informational alternatives and are therefore
    never individually required.
    """
    statuses: list[ToolStatus] = []
    for tool in _TOOL_GROUPS:
        interfaces = [_interface_status(tool, interface) for interface in tool.interfaces]
        usable = [status for status in interfaces if status.available]
        version = next((status.version for status in usable if status.version), None)
        summary = ToolStatus(
            name=tool.name,
            kind="tool",
            required=tool.required,
            available=bool(usable),
            detail=(
                "usable via " + ", ".join(_interface_label(tool, status) for status in usable)
                if usable
                else tool.install_hint
            ),
            version=version,
        )
        statuses.extend([summary, *interfaces])
    return statuses


def _interface_label(tool: _ToolGroup, status: ToolStatus) -> str:
    prefix = f"{tool.name} ("
    if status.name.startswith(prefix) and status.name.endswith(")"):
        return status.name[len(prefix) : -1]
    return status.kind


def _interface_status(tool: _ToolGroup, interface: _InterfaceCheck) -> ToolStatus:
    if interface.kind == "python-module":
        available, version = python_import_check(interface.probe_names[0])
        detail = "importable" if available else "not importable"
    else:
        found = next(
            (
                (candidate, path)
                for candidate in interface.probe_names
                if (path := which_executable(candidate)) is not None
            ),
            None,
        )
        available = found is not None
        if found is not None:
            candidate, path = found
            detail = path
            version = executable_version(candidate, args=list(interface.version_args))
        else:
            detail = f"not on PATH (tried: {', '.join(interface.probe_names)})"
            version = None
    return ToolStatus(
        name=f"{tool.name} ({interface.label})",
        kind=interface.kind,
        required=False,
        available=available,
        detail=detail,
        version=version,
        group=tool.name,
    )


def _executable_statuses() -> list[ToolStatus]:
    statuses: list[ToolStatus] = []
    for executable, config in _EXECUTABLE_CHECKS.items():
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
