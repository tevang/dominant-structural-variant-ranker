from __future__ import annotations

import inspect
import shutil
import tempfile
from pathlib import Path
from typing import Any

from rdkit import Chem

from dsvr.runners.subprocess_utils import run_command


class MolscrubUnavailableError(RuntimeError):
    """Raised when neither molscrub Python API nor CLI is available."""


class MolscrubExecutionError(RuntimeError):
    """Raised when molscrub is available but candidate generation fails."""


def inspect_molscrub() -> dict[str, Any]:
    api_available = False
    api_signature = None
    api_error = None
    try:
        from molscrub import Scrub  # type: ignore[import-not-found]

        api_available = True
        api_signature = str(inspect.signature(Scrub))
    except Exception as exc:  # pragma: no cover - environment dependent
        api_error = f"{type(exc).__name__}: {exc}"
    return {
        "python_api_available": api_available,
        "python_api_signature": api_signature,
        "python_api_error": api_error,
        "scrub_py": shutil.which("scrub.py"),
        "molscrub_executable": shutil.which("molscrub"),
    }


def generate_molscrub_candidates(
    molecule: Chem.Mol,
    *,
    ph_low: float,
    ph_high: float,
) -> tuple[list[Chem.Mol], str, str]:
    try:
        return generate_molscrub_candidates_python(molecule, ph_low=ph_low, ph_high=ph_high)
    except MolscrubUnavailableError:
        return generate_molscrub_candidates_cli(molecule, ph_low=ph_low, ph_high=ph_high)


def generate_molscrub_candidates_python(
    molecule: Chem.Mol,
    *,
    ph_low: float,
    ph_high: float,
) -> tuple[list[Chem.Mol], str, str]:
    try:
        from molscrub import Scrub  # type: ignore[import-not-found]
    except Exception as exc:
        raise MolscrubUnavailableError(_install_message()) from exc

    kwargs: dict[str, Any] = {"ph_low": ph_low, "ph_high": ph_high}
    signature = inspect.signature(Scrub)
    if "skip_gen3d" in signature.parameters:
        kwargs["skip_gen3d"] = True
    if "skip_gen3d" not in signature.parameters and "gen3d" in signature.parameters:
        kwargs["gen3d"] = False

    try:
        scrubber = Scrub(**kwargs)
        candidates = list(scrubber(molecule))
    except Exception as exc:
        raise MolscrubExecutionError(f"molscrub Python API failed: {exc}") from exc

    return candidates, "molscrub-python-api", f"Scrub({kwargs!r})(mol)"


def generate_molscrub_candidates_cli(
    molecule: Chem.Mol,
    *,
    ph_low: float,
    ph_high: float,
) -> tuple[list[Chem.Mol], str, str]:
    executable = shutil.which("scrub.py") or shutil.which("molscrub")
    if executable is None:
        raise MolscrubUnavailableError(_install_message())

    with tempfile.TemporaryDirectory(prefix="dsvr_molscrub_") as tmpdir_raw:
        tmpdir = Path(tmpdir_raw)
        input_sdf = tmpdir / "input.sdf"
        output_sdf = tmpdir / "output.sdf"
        writer = Chem.SDWriter(str(input_sdf))
        writer.write(molecule)
        writer.close()

        command = [
            executable,
            str(input_sdf),
            "-o",
            str(output_sdf),
            "--ph",
            str(ph_low if ph_low == ph_high else (ph_low + ph_high) / 2.0),
        ]
        if _cli_help_mentions(executable, "--skip_gen3d"):
            command.append("--skip_gen3d")
        completed = run_command(
            command,
            timeout_s=300,
            command_name="molscrub",
            check=False,
        )
        if completed.returncode != 0:
            raise MolscrubExecutionError(
                "molscrub CLI failed with exit code "
                f"{completed.returncode}: {completed.stdout.strip()}"
            )
        if not output_sdf.exists():
            raise MolscrubExecutionError("molscrub CLI did not create expected output SDF")
        supplier = Chem.SDMolSupplier(str(output_sdf), sanitize=True, removeHs=False)
        candidates = [mol for mol in supplier if mol is not None]
        return candidates, "molscrub-cli", " ".join(command)


def _cli_help_mentions(executable: str, option: str) -> bool:
    completed = run_command(
        [executable, "-h"],
        timeout_s=15,
        command_name="molscrub_help",
        check=False,
    )
    return option in completed.stdout


def _install_message() -> str:
    return (
        "molscrub is required for protomer candidate generation but is not installed. "
        "Install it in the active environment, for example with "
        "`pip install git+https://github.com/forlilab/molscrub.git`, then run "
        "`dsvr doctor` to verify availability."
    )
