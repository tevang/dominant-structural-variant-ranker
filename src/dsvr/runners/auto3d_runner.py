from __future__ import annotations

import importlib.metadata
import importlib.util
import logging
import shutil
import subprocess
import sys
import uuid
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem

from dsvr.runners.subprocess_utils import ExternalToolError, run_command

_LOGGER = logging.getLogger("dsvr.auto3d")

# Error classes shared by stage-level failure bookkeeping (auto3d-failure-reporting).
CLASS_ENGINE_INCOMPATIBLE = "ENGINE_INCOMPATIBLE"
CLASS_INFRA_MULTIPROCESSING = "INFRA_MULTIPROCESSING"
CLASS_CUDA_UNAVAILABLE = "CUDA_UNAVAILABLE"
CLASS_TIMEOUT = "TIMEOUT"
CLASS_EXECUTION = "EXECUTION_ERROR"

def classify_auto3d_failure(text: str) -> str:
    """Classify an Auto3D failure text into a stable error class."""

    if any(
        f"Only {engine} can handle" in text
        for engine in ("AIMNET", "AIMNet2", "ANI2x", "ANI2xt")
    ):
        # Auto3D validation message, e.g. "Only AIMNET can handle: [...]"
        return CLASS_ENGINE_INCOMPATIBLE
    if "SemLock" in text and "fork context" in text:
        return CLASS_INFRA_MULTIPROCESSING
    cuda_signatures = (
        "no CUDA-capable device",
        "CUDA initialization",
        "CUDA driver",
        "Found no NVIDIA driver",
        "cuInit",
    )
    if any(signature in text for signature in cuda_signatures):
        return CLASS_CUDA_UNAVAILABLE
    if "timed out after" in text:
        return CLASS_TIMEOUT
    return CLASS_EXECUTION


class Auto3DUnavailableError(RuntimeError):
    """Raised when Auto3D is required but unavailable."""


class Auto3DExecutionError(RuntimeError):
    """Raised when Auto3D execution fails."""


@dataclass(frozen=True)
class _Auto3DCachePaths:
    xdg_cache_home: Path
    warp_cache_path: Path
    aimnet_cache_dir: Path


def inspect_auto3d() -> dict[str, str | bool | None]:
    return {
        "python_api_available": importlib.util.find_spec("Auto3D") is not None,
        "auto3d_executable": shutil.which("auto3d"),
        "auto3D_executable": shutil.which("auto3D"),
        "Auto3D_executable": shutil.which("Auto3D"),
    }


# --------------------------------------------------------------------------
# Engine capability awareness
#
# Static engine→capability table derived from Auto3D's own validation:
# ANI2x/ANI2xt handle neutral, closed-shell molecules over their training
# elements; AIMNET/AIMNet2 accept a broader main-group set plus charged and
# open-shell species ("Only AIMNET can handle" otherwise). The table can be
# overridden via ``config.auto3d.engine_element_overrides`` if Auto3D's real
# validation drifts from it.
# --------------------------------------------------------------------------

_AUTO3D_ENGINE_ELEMENTS: dict[str, frozenset[str]] = {
    "ANI2x": frozenset({"H", "C", "N", "O", "F", "S", "Cl"}),
    "ANI2xt": frozenset({"H", "C", "N", "O", "F", "S", "Cl"}),
    "AIMNET": frozenset({"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "As", "Se", "Br", "I"}),
    "AIMNet2": frozenset({"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "As", "Se", "Br", "I"}),
}

_ENGINES_ALLOWING_NONNEUTRAL: frozenset[str] = frozenset({"AIMNET", "AIMNet2", "auto"})


def engine_supports_molecule(
    engine: str,
    mol: Chem.Mol,
    *,
    element_overrides: dict[str, list[str]] | None = None,
) -> bool:
    """Return True when ``engine`` can handle ``mol``'s composition.

    Unknown engines are assumed capable (Auto3D's own validation remains the
    fallback arbiter). ``auto`` routes internally inside Auto3D.
    """

    if engine == "auto":
        return True
    allowed: frozenset[str] | None
    if element_overrides and engine in element_overrides:
        allowed = frozenset(element_overrides[engine])
    else:
        allowed = _AUTO3D_ENGINE_ELEMENTS.get(engine)
    if allowed is None:
        return True
    symbols = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if not symbols <= allowed:
        return False
    if engine not in _ENGINES_ALLOWING_NONNEUTRAL:
        if Chem.GetFormalCharge(mol) != 0:
            return False
        if any(atom.GetNumRadicalElectrons() > 0 for atom in mol.GetAtoms()):
            return False
    return True


def partition_by_engine(
    items: list[tuple[Any, Chem.Mol]],
    engines: list[str],
    *,
    element_overrides: dict[str, list[str]] | None = None,
) -> tuple[dict[str, list[Any]], list[Any]]:
    """Route every item to the first configured engine that supports its molecule.

    ``items`` are ``(payload, mol)`` pairs (payload is typically a record id).
    ``engines`` is the ordered preference list. Returns ``(assignments,
    incompatible)``: assignments maps each engine to its payloads in input
    order; incompatible lists payloads no configured engine supports — those
    must go straight to the recorded RDKit fallback instead of being offered
    to Auto3D repeatedly.
    """

    assignments: dict[str, list[Any]] = {}
    incompatible: list[Any] = []
    for payload, mol in items:
        engine = next(
            (
                candidate
                for candidate in engines
                if engine_supports_molecule(candidate, mol, element_overrides=element_overrides)
            ),
            None,
        )
        if engine is None:
            incompatible.append(payload)
        else:
            assignments.setdefault(engine, []).append(payload)
    return assignments, incompatible


def run_auto3d(
    input_path: Path,
    output_dir: Path,
    *,
    k: int,
    model: str,
    internal_tautomer_stereo_enum: bool,
    mpi_np: int | None = None,
    cpu_workers: int | None = None,
    memory_gb: int | None = None,
    capacity: int | None = None,
    max_confs: int | None = None,
    patience: int | None = None,
    threshold: float | None = None,
    opt_steps: int | None = None,
    use_gpu: bool = False,
    stream_output: bool = False,
    timeout_s: int | None = None,
    isomer_enum_only: bool = False,
) -> tuple[Path, list[str]]:
    executable = _find_executable()
    python_api_available = importlib.util.find_spec("Auto3D") is not None
    if executable is None and not python_api_available:
        raise Auto3DUnavailableError(
            "Auto3D is required for Auto3D seeding but is not installed. Install it with "
            "`pip install Auto3D` or `conda install -c conda-forge auto3d`, then run "
            "`dsvr doctor` to verify availability."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    cache_paths = _auto3d_cache_paths()
    auto3d_major = _auto3d_major_version()
    v3_wrapper_script = (
        _ensure_auto3d_v3_wrapper_script(output_dir)
        if auto3d_major is not None and auto3d_major >= 3
        else None
    )
    wrapper_script = (
        _ensure_auto3d_wrapper_script(output_dir)
        if auto3d_major is not None and auto3d_major < 3
        else None
    )
    output_sdf = output_dir / "auto3d_output.sdf"
    job_name_base = f"{_output_dir_name(output_sdf)}_{uuid.uuid4().hex[:8]}"
    failures: list[str] = []
    use_gpu = _should_use_gpu(use_gpu)
    for command in _command_candidates(
        executable,
        v3_wrapper_script,
        wrapper_script,
        input_path,
        output_sdf,
        job_name_base,
        k=k,
        model=model,
        internal_tautomer_stereo_enum=internal_tautomer_stereo_enum,
        isomer_enum_only=isomer_enum_only,
        mpi_np=mpi_np,
        cpu_workers=cpu_workers,
        memory_gb=memory_gb,
        capacity=capacity,
        max_confs=max_confs,
        patience=patience,
        threshold=threshold,
        opt_steps=opt_steps,
        use_gpu=use_gpu,
    ):
        try:
            completed = run_command(
                command,
                cwd=output_dir,
                timeout_s=timeout_s,
                log_dir=output_dir / "logs",
                command_name="auto3d",
                env=_auto3d_env(
                    cache_paths,
                    mpi_np=mpi_np,
                    cpu_workers=cpu_workers,
                    use_gpu=use_gpu,
                ),
                stream_output=stream_output,
                check=False,
            )
        except ExternalToolError as exc:
            if not (exc.metadata or {}).get("timed_out"):
                raise
            # A timeout is a property of the workload and the limit, not of
            # the command form — retrying every command candidate would burn
            # the same limit again. Fail the invocation once; the caller's
            # smaller-batch retry path decides what to do next.
            raise Auto3DExecutionError(
                f"Auto3D timed out after {timeout_s} s: {' '.join(command)}"
            ) from exc
        if completed.returncode != 0:
            guessed = _find_output_sdf(
                output_dir, input_path=input_path, job_name=job_name_base
            )
            if guessed is not None and _sdf_contains_records(guessed):
                # Auto3D exits nonzero when some inputs produce no output
                # (e.g. charged molecules it cannot optimize), but it still
                # writes usable results for the rest. Treat that as a partial
                # success: downstream stages fill the missing variants with
                # the RDKit fallback and record per-variant warnings.
                return guessed, command
            failure_text = _completed_output_tail(completed)
            failures.append(
                f"{' '.join(command)} exited {completed.returncode}: "
                f"{failure_text}"
            )
            if _is_terminal_auto3d_selection_failure(failure_text):
                break
            continue
        if output_sdf.exists():
            return output_sdf, command
        guessed = _find_output_sdf(
            output_dir, input_path=input_path, job_name=job_name_base
        )
        if guessed is not None:
            return guessed, command
        failures.append(
            f"{' '.join(command)} exited 0 but did not produce an output SDF"
        )

    raise Auto3DExecutionError("Auto3D failed. Tried commands:\n" + "\n".join(failures))


def _completed_output_tail(
    completed: object,
    *,
    max_chars: int = 4000,
) -> str:
    stdout = getattr(completed, "stdout", "") or ""
    stderr = getattr(completed, "stderr", "") or ""
    text = (str(stdout) + "\n" + str(stderr)).strip()
    if len(text) <= max_chars:
        return text
    return text[-max_chars:]


def _is_terminal_auto3d_selection_failure(text: str) -> bool:
    """Return True when retrying a legacy Auto3D invocation is unlikely to help."""

    # Engine-incompatibility validation ("Only AIMNET can handle: ...") and
    # CUDA-absence errors cannot disappear across command forms; the caller
    # must switch engine or mode instead of retrying.
    error_class = classify_auto3d_failure(text)
    if error_class in {CLASS_ENGINE_INCOMPATIBLE, CLASS_CUDA_UNAVAILABLE}:
        return True
    if "Dropped(Oscillating)" in text and "Converged: 0" in text:
        return True
    return "reorder_sdf" in text and "Invalid input file" in text and "_out.sdf" in text


def _find_executable() -> str | None:
    for executable in ("auto3d", "auto3D", "Auto3D"):
        path = shutil.which(executable)
        if path is not None:
            return path
    return None


def _command_candidates(
    executable: str | None,
    v3_wrapper_script: Path | None,
    wrapper_script: Path | None,
    input_path: Path,
    output_sdf: Path,
    job_name_base: str,
    *,
    k: int,
    model: str,
    internal_tautomer_stereo_enum: bool,
    isomer_enum_only: bool,
    mpi_np: int | None,
    cpu_workers: int | None,
    memory_gb: int | None,
    capacity: int | None,
    max_confs: int | None,
    patience: int | None,
    threshold: float | None,
    opt_steps: int | None,
    use_gpu: bool,
) -> list[list[str]]:
    commands: list[list[str]] = []
    auto3d_major = _auto3d_major_version()
    if auto3d_major is not None and auto3d_major >= 3 and (
        v3_wrapper_script is not None or executable is not None
    ):
        v3_prefix: list[str]
        if v3_wrapper_script is not None:
            v3_prefix = [sys.executable, str(v3_wrapper_script.resolve())]
        else:
            v3_prefix = [str(executable)]

        v3 = [
            *v3_prefix,
            "run",
            str(input_path.resolve()),
            "--k",
            str(k),
            "--job-name",
            f"{job_name_base}_v3",
            "--engine",
            _auto3d_v3_engine(model),
            "--tauto-engine",
            "rdkit",
            "--isomer-engine",
            "rdkit",
        ]
        if isomer_enum_only:
            # Isomer enumeration without tautomer re-enumeration: used by the
            # auto3d_enumerate stereo policy for unspecified-stereo sub-batches.
            v3.extend(["--no-enumerate-tautomer", "--enumerate-isomer"])
        elif internal_tautomer_stereo_enum:
            v3.extend(["--enumerate-tautomer", "--enumerate-isomer"])
        else:
            v3.extend(["--no-enumerate-tautomer", "--no-enumerate-isomer"])
        v3.append("--gpu" if use_gpu else "--no-gpu")
        commands.append(v3)

    if wrapper_script is not None:
        wrapper = [
            sys.executable,
            str(wrapper_script.resolve()),
            str(input_path.resolve()),
            "--k",
            str(k),
            "--job_name",
            f"{job_name_base}_shim",
            "--optimizing_engine",
            model,
            "--isomer_engine",
            "rdkit",
            "--tauto_engine",
            "rdkit",
            "--enumerate_tautomer",
            "False" if isomer_enum_only else ("True" if internal_tautomer_stereo_enum else "False"),
            "--enumerate_isomer",
            "True" if (internal_tautomer_stereo_enum or isomer_enum_only) else "False",
        ]
        if mpi_np is not None:
            wrapper.extend(["--mpi_np", str(mpi_np)])
        if memory_gb is not None:
            wrapper.extend(["--memory", str(memory_gb)])
        if capacity is not None:
            wrapper.extend(["--capacity", str(capacity)])
        commands.append(wrapper)

    if executable is not None and (auto3d_major is None or auto3d_major < 3):
        commands.extend(
            [
                [
                    executable,
                    str(input_path.resolve()),
                    "--k",
                    str(k),
                    "--job_name",
                    f"{job_name_base}_legacy1",
                    "--optimizing_engine",
                    model,
                    "--isomer_engine",
                    "rdkit",
                    "--tauto_engine",
                    "rdkit",
                ],
                [
                    executable,
                    str(input_path.resolve()),
                    "--k",
                    str(k),
                    "--job_name",
                    f"{job_name_base}_legacy2",
                    "--enumerate_tautomer",
                    "False" if isomer_enum_only else ("True" if internal_tautomer_stereo_enum else "False"),
                    "--enumerate_isomer",
                    "True" if (internal_tautomer_stereo_enum or isomer_enum_only) else "False",
                    "--optimizing_engine",
                    model,
                    "--isomer_engine",
                    "rdkit",
                    "--tauto_engine",
                    "rdkit",
                ],
            ]
        )
        for command in commands:
            if _is_auto3d_v3_run(command):
                continue
            if command and command[0] == sys.executable:
                continue
            if mpi_np is not None:
                command.extend(["--mpi_np", str(mpi_np)])
            if memory_gb is not None:
                command.extend(["--memory", str(memory_gb)])
            if capacity is not None:
                command.extend(["--capacity", str(capacity)])
    if max_confs is not None:
        for command in commands:
            if _is_auto3d_v3_run(command):
                command.extend(["--max-confs", str(max_confs)])
            else:
                command.extend(["--max_confs", str(max_confs)])
    if patience is not None:
        for command in commands:
            command.extend(["--patience", str(patience)])
    if threshold is not None:
        for command in commands:
            command.extend(["--threshold", str(threshold)])
    if opt_steps is not None:
        for command in commands:
            if _is_auto3d_v3_run(command):
                command.extend(["--opt-steps", str(opt_steps)])
            else:
                command.extend(["--opt_steps", str(opt_steps)])
    for command in commands:
        if not _is_auto3d_v3_run(command):
            command.extend(["--use_gpu", "True" if use_gpu else "False"])
    return commands


def _auto3d_env(
    cache_paths: _Auto3DCachePaths,
    *,
    mpi_np: int | None,
    cpu_workers: int | None,
    use_gpu: bool,
) -> dict[str, str]:
    """Environment overrides to keep Auto3D/AIMNET/Warp caches repo-local.

    CPU mode hides CUDA devices (``CUDA_VISIBLE_DEVICES=""``). Warp may still
    emit a one-line CUDA-related notice at import even with devices hidden —
    recent Warp versions treat an empty device list as a clean CPU fallback,
    so the residual line is informational and not an initialization error.
    No GPU-related command-line flags are passed in CPU mode.
    """

    env: dict[str, str] = {
        "XDG_CACHE_HOME": str(cache_paths.xdg_cache_home),
        "WARP_CACHE_PATH": str(cache_paths.warp_cache_path),
        "AIMNET_CACHE_DIR": str(cache_paths.aimnet_cache_dir),
    }
    if not use_gpu:
        env["CUDA_VISIBLE_DEVICES"] = ""

    if mpi_np is None:
        return env

    per_process_threads = mpi_np
    if cpu_workers is not None and cpu_workers > 1:
        per_process_threads = max(1, mpi_np // cpu_workers)
    value = str(per_process_threads)
    env.update(
        {
            "OMP_NUM_THREADS": value,
            "MKL_NUM_THREADS": value,
            "OPENBLAS_NUM_THREADS": value,
            "NUMEXPR_NUM_THREADS": value,
        }
    )
    return env


def _cpu_worker_indices(cpu_workers: int) -> str:
    return ",".join(str(index) for index in range(cpu_workers))


def _find_output_sdf(
    output_dir: Path,
    *,
    input_path: Path | None = None,
    job_name: str | None = None,
) -> Path | None:
    candidates: list[Path] = []
    if job_name:
        # Prefer outputs created by this invocation: its unique job name
        # appears in the job directory name, so a scoped search cannot pick
        # up a stale ``*_out.sdf`` left behind by an earlier crashed run in
        # the same persistent output_dir (e.g. the final-3D stage).
        candidates.extend(sorted(output_dir.glob(f"*{job_name}*/**/*_out.sdf")))
        candidates.extend(sorted(output_dir.glob(f"*{job_name}*/**/*_3d.sdf")))
    if not candidates:
        candidates = sorted(
            list(output_dir.glob("**/*_out.sdf")) + list(output_dir.glob("**/*_3d.sdf"))
        )
    if input_path is not None and job_name:
        # Auto3D creates job directories next to the input file
        # (<input_stem>_<job_name>/...), which can live outside output_dir
        # when callers pass a sibling directory as output_dir (e.g. the
        # tautomer-filtering stage). Search there too, restricted to this
        # invocation's unique job name to avoid matching stale outputs.
        matches = list(input_path.parent.glob(f"*{job_name}*/**/*_out.sdf"))
        matches += input_path.parent.glob(f"*{job_name}*/**/*_3d.sdf")
        candidates.extend(sorted(matches))
    return candidates[0] if candidates else None


def _sdf_contains_records(path: Path) -> bool:
    try:
        return b"$$$$" in path.read_bytes()
    except OSError:
        return False


def _python_executable() -> str:
    return shutil.which("python") or shutil.which("python3") or "python"


def _output_dir_name(output_sdf: Path) -> str:
    return output_sdf.parent.name


def _repo_root() -> Path:
    # auto3d_runner.py -> runners -> dsvr -> src -> repo root
    return Path(__file__).resolve().parents[3]


def _auto3d_cache_paths() -> _Auto3DCachePaths:
    root = _repo_root()
    uv_root = root / ".uv"
    paths = _Auto3DCachePaths(
        xdg_cache_home=uv_root / "xdg-cache",
        warp_cache_path=uv_root / "warp-cache",
        aimnet_cache_dir=uv_root / "aimnet-cache",
    )
    paths.xdg_cache_home.mkdir(parents=True, exist_ok=True)
    paths.warp_cache_path.mkdir(parents=True, exist_ok=True)
    paths.aimnet_cache_dir.mkdir(parents=True, exist_ok=True)
    return paths


_AUTO3D_WARNING_OBSERVERS: list[Callable[[str], None]] = []
_GPU_PROBE_CACHE: dict[str, bool] = {}
_EMITTED_ONCE_KEYS: set[str] = set()


def register_auto3d_warning_observer(observer: Callable[[str], None]) -> None:
    """Register a sink receiving one-time Auto3D environment notices.

    The workflow engine registers ``ProgressRecorder.warning`` here so that
    environment degradations (e.g. GPU unusable) land in ``warnings.jsonl``
    exactly once per run instead of once per molecule.
    """

    _AUTO3D_WARNING_OBSERVERS.append(observer)


def _emit_auto3d_notice_once(key: str, message: str) -> None:
    if key in _EMITTED_ONCE_KEYS:
        return
    _EMITTED_ONCE_KEYS.add(key)
    _LOGGER.warning(message)
    for observer in _AUTO3D_WARNING_OBSERVERS:
        observer(message)


def _should_use_gpu(requested: bool) -> bool:
    if not requested:
        return False
    if _probe_gpu_usable():
        return True
    _emit_auto3d_notice_once(
        "gpu-unusable",
        "Auto3D GPU use was requested, but the runtime GPU probe failed "
        f"({_GPU_PROBE_CACHE.get('reason', 'no usable CUDA device')}); "
        "degrading to CPU mode for all Auto3D stages of this run.",
    )
    return False


def _probe_gpu_usable() -> bool:
    """Runtime GPU usability verdict, cached per process (one probe per run).

    Replaces the previous ``/dev/nvidia*`` device-node check: device nodes
    can be present while the driver/userspace stack is broken, so usability
    is verified by actually initializing CUDA in a probe subprocess.
    """

    if "verdict" in _GPU_PROBE_CACHE:
        return _GPU_PROBE_CACHE["verdict"]

    verdict = False
    reason = "no usable CUDA device"
    smi = shutil.which("nvidia-smi")
    smi_ok = False
    if smi is not None:
        try:
            smi_ok = (
                subprocess.run(
                    [smi],
                    check=False,
                    capture_output=True,
                    timeout=30,
                ).returncode
                == 0
            )
        except (OSError, subprocess.TimeoutExpired):
            smi_ok = False
    if smi is None or not smi_ok:
        reason = "nvidia-smi missing or failed; driver stack unusable"
    else:
        try:
            probe = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    "import sys, torch; sys.exit(0 if torch.cuda.is_available() else 3)",
                ],
                check=False,
                capture_output=True,
                timeout=180,
            )
            verdict = probe.returncode == 0
            if not verdict:
                reason = "torch.cuda.is_available() is False in the Auto3D interpreter"
        except (OSError, subprocess.TimeoutExpired) as exc:
            reason = f"CUDA probe subprocess failed: {exc}"

    _GPU_PROBE_CACHE["verdict"] = verdict
    _GPU_PROBE_CACHE["reason"] = "" if verdict else reason
    return verdict


def _auto3d_major_version() -> int | None:
    try:
        raw = importlib.metadata.version("auto3d")
    except importlib.metadata.PackageNotFoundError:
        return None
    try:
        return int(raw.split(".", maxsplit=1)[0])
    except ValueError:
        return None


def _is_auto3d_v3_run(command: list[str]) -> bool:
    if not command:
        return False
    # Support both direct invocation:
    #   auto3d run <input> ...
    # and our sandbox-tolerant wrapper:
    #   python _auto3d_v3_wrapper.py run <input> ...
    return (len(command) > 1 and command[1] == "run") or (
        len(command) > 2 and command[2] == "run"
    )


def _auto3d_v3_engine(model: str) -> str:
    normalized = model.strip()
    mapping = {
        "AIMNet2": "aimnet2",
        "AIMNET2": "aimnet2",
        "AIMNET": "AIMNET",
        "ANI2x": "ANI2x",
        "ANI2xt": "ANI2xt",
        "auto": "auto",
        "AUTO": "auto",
    }
    return mapping.get(normalized, normalized)


def _ensure_auto3d_v3_wrapper_script(output_dir: Path) -> Path:
    """Wrap Auto3D v3 CLI to tolerate restricted socket syscalls.

    Some sandboxed environments block ``socket.setsockopt`` used by
    ``multiprocessing.managers.SyncManager``. Auto3D v3 can trigger Manager
    creation even for simple CPU runs, so we patch ``setsockopt`` to ignore
    PermissionError and keep the CLI usable.
    """

    script = output_dir / "_auto3d_v3_wrapper.py"
    script.write_text(
        """
from __future__ import annotations

import socket


_orig_setsockopt = socket.socket.setsockopt


def _patched_setsockopt(self, *args, **kwargs):
    try:
        return _orig_setsockopt(self, *args, **kwargs)
    except PermissionError:
        return None


socket.socket.setsockopt = _patched_setsockopt

from dsvr.runners._auto3d_mp_shim import install_fake_manager

install_fake_manager()

from Auto3D.presentation.auto3Dcli import cli


if __name__ == "__main__":
    raise SystemExit(cli())
""".lstrip(),
        encoding="utf-8",
    )
    return script


def _ensure_auto3d_wrapper_script(output_dir: Path) -> Path:
    script = output_dir / "_auto3d_wrapper.py"
    script.write_text(
        """
from __future__ import annotations

import importlib.metadata as _md
import sys
import types


class _DistributionNotFound(Exception):
    pass


def _get_distribution(name: str):
    try:
        return types.SimpleNamespace(version=_md.version(name))
    except _md.PackageNotFoundError as exc:
        raise _DistributionNotFound(name) from exc


_pkg_resources = types.ModuleType("pkg_resources")
_pkg_resources.get_distribution = _get_distribution
_pkg_resources.DistributionNotFound = _DistributionNotFound
sys.modules["pkg_resources"] = _pkg_resources

from dsvr.runners._auto3d_mp_shim import install_fake_manager

install_fake_manager()

import Auto3D.auto3D as _auto3d_module


_original_isomer_wrapper = _auto3d_module.isomer_wraper


def _nonempty_chunk_info(chunk_info):
    filtered = []
    for path, workdir in chunk_info:
        try:
            has_records = any(line.strip() for line in open(path, encoding="utf-8"))
        except OSError:
            has_records = True
        if has_records:
            filtered.append((path, workdir))
    return filtered


def _isomer_wrapper_skip_empty_chunks(chunk_info, args, queue, logging_queue):
    filtered = _nonempty_chunk_info(chunk_info)
    if filtered:
        return _original_isomer_wrapper(filtered, args, queue, logging_queue)
    done_count = 1 if isinstance(args.gpu_idx, int) else len(args.gpu_idx)
    for _ in range(done_count):
        queue.put("Done")
    return None


_auto3d_module.isomer_wraper = _isomer_wrapper_skip_empty_chunks

from Auto3D.auto3Dcli import cli


if __name__ == "__main__":
    cli()
""".lstrip(),
        encoding="utf-8",
    )
    return script
