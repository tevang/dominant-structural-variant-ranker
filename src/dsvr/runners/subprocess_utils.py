from __future__ import annotations

import subprocess
from pathlib import Path


def run_command(
    command: list[str],
    cwd: Path | None = None,
    timeout_s: int | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        timeout=timeout_s,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
