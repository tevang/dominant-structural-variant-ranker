from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

from dsvr.runners.subprocess_utils import run_command


class Auto3DUnavailableError(RuntimeError):
    """Raised when Auto3D is required but unavailable."""


class Auto3DExecutionError(RuntimeError):
    """Raised when Auto3D execution fails."""


def inspect_auto3d() -> dict[str, str | bool | None]:
    return {
        "python_api_available": importlib.util.find_spec("Auto3D") is not None,
        "auto3d_executable": shutil.which("auto3d"),
        "auto3D_executable": shutil.which("auto3D"),
        "Auto3D_executable": shutil.which("Auto3D"),
    }


def run_auto3d(
    input_path: Path,
    output_dir: Path,
    *,
    k: int,
    model: str,
    internal_tautomer_stereo_enum: bool,
) -> tuple[Path, list[str]]:
    executable = _find_executable()
    if executable is None:
        raise Auto3DUnavailableError(
            "Auto3D is required for Auto3D seeding but is not installed. Install it with "
            "`pip install Auto3D` or `conda install -c conda-forge auto3d`, then run "
            "`dsvr doctor` to verify availability."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    output_sdf = output_dir / "auto3d_output.sdf"
    failures: list[str] = []
    for command in _command_candidates(
        executable,
        input_path,
        output_sdf,
        k=k,
        model=model,
        internal_tautomer_stereo_enum=internal_tautomer_stereo_enum,
    ):
        completed = run_command(
            command,
            cwd=output_dir,
            timeout_s=None,
            log_dir=output_dir / "logs",
            command_name="auto3d",
            check=False,
        )
        if completed.returncode != 0:
            failures.append(
                f"{' '.join(command)} exited {completed.returncode}: "
                f"{completed.stdout.strip()}"
            )
            continue
        if output_sdf.exists():
            return output_sdf, command
        guessed = _find_output_sdf(output_dir)
        if guessed is not None:
            return guessed, command
        failures.append(
            f"{' '.join(command)} exited 0 but did not produce an output SDF"
        )

    raise Auto3DExecutionError("Auto3D failed. Tried commands:\n" + "\n".join(failures))


def _find_executable() -> str | None:
    for executable in ("auto3d", "auto3D", "Auto3D"):
        path = shutil.which(executable)
        if path is not None:
            return path
    return None


def _command_candidates(
    executable: str,
    input_path: Path,
    output_sdf: Path,
    *,
    k: int,
    model: str,
    internal_tautomer_stereo_enum: bool,
) -> list[list[str]]:
    commands = [
        [
            executable,
            "run",
            str(input_path),
            "--k",
            str(k),
            "--engine",
            model,
            "--out",
            str(output_sdf),
        ],
        [
            executable,
            "--path",
            str(input_path),
            "--k",
            str(k),
            "--optimizing_engine",
            model,
            "--output",
            str(output_sdf),
        ],
    ]
    if not internal_tautomer_stereo_enum:
        commands[0].extend(["--no-tautomer", "--no-isomer"])
        commands[1].extend(
            [
                "--enumerate_tautomer",
                "False",
                "--enumerate_isomer",
                "False",
            ]
        )
    return commands


def _find_output_sdf(output_dir: Path) -> Path | None:
    candidates = sorted(output_dir.glob("*.sdf"))
    return candidates[0] if candidates else None
