from __future__ import annotations

import importlib.util
import shutil

from dsvr.models import ToolStatus

PYTHON_MODULES = ("rdkit", "molscrub", "Auto3D", "psi4", "pyscf")
EXECUTABLES = ("xtb", "crest", "censo", "psi4")


def check_tools() -> list[ToolStatus]:
    statuses: list[ToolStatus] = []
    for module_name in PYTHON_MODULES:
        available = importlib.util.find_spec(module_name) is not None
        statuses.append(
            ToolStatus(
                name=module_name,
                kind="python-module",
                required=False,
                available=available,
                detail="importable" if available else "not importable",
            )
        )
    for executable in EXECUTABLES:
        path = shutil.which(executable)
        statuses.append(
            ToolStatus(
                name=executable,
                kind="executable",
                required=False,
                available=path is not None,
                detail=path or "not on PATH",
            )
        )
    return statuses

