from __future__ import annotations

import importlib.util
import shutil
import uuid
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
) -> tuple[Path, list[str]]:
    executable = _find_executable()
    if executable is None:
        raise Auto3DUnavailableError(
            "Auto3D is required for Auto3D seeding but is not installed. Install it with "
            "`pip install Auto3D` or `conda install -c conda-forge auto3d`, then run "
            "`dsvr doctor` to verify availability."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    wrapper_script = _ensure_auto3d_wrapper_script(output_dir)
    output_sdf = output_dir / "auto3d_output.sdf"
    job_name_base = f"{_output_dir_name(output_sdf)}_{uuid.uuid4().hex[:8]}"
    failures: list[str] = []
    for command in _command_candidates(
        executable,
        wrapper_script,
        input_path,
        output_sdf,
        job_name_base,
        k=k,
        model=model,
        internal_tautomer_stereo_enum=internal_tautomer_stereo_enum,
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
        completed = run_command(
            command,
            cwd=output_dir,
            timeout_s=None,
            log_dir=output_dir / "logs",
            command_name="auto3d",
            env=_auto3d_env(mpi_np, cpu_workers),
            stream_output=stream_output,
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
    wrapper_script: Path,
    input_path: Path,
    output_sdf: Path,
    job_name_base: str,
    *,
    k: int,
    model: str,
    internal_tautomer_stereo_enum: bool,
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
    wrapper = [
        _python_executable(),
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
        "True" if internal_tautomer_stereo_enum else "False",
        "--enumerate_isomer",
        "True" if internal_tautomer_stereo_enum else "False",
    ]
    if mpi_np is not None:
        wrapper.extend(["--mpi_np", str(mpi_np)])
    if cpu_workers is not None and cpu_workers > 1:
        wrapper.extend(["--gpu_idx", _cpu_worker_indices(cpu_workers)])
    if memory_gb is not None:
        wrapper.extend(["--memory", str(memory_gb)])
    if capacity is not None:
        wrapper.extend(["--capacity", str(capacity)])
    commands = [
        wrapper,
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
            "True" if internal_tautomer_stereo_enum else "False",
            "--enumerate_isomer",
            "True" if internal_tautomer_stereo_enum else "False",
            "--optimizing_engine",
            model,
            "--isomer_engine",
            "rdkit",
            "--tauto_engine",
            "rdkit",
        ],
    ]
    for command in commands[1:]:
        if mpi_np is not None:
            command.extend(["--mpi_np", str(mpi_np)])
        if cpu_workers is not None and cpu_workers > 1:
            command.extend(["--gpu_idx", _cpu_worker_indices(cpu_workers)])
        if memory_gb is not None:
            command.extend(["--memory", str(memory_gb)])
        if capacity is not None:
            command.extend(["--capacity", str(capacity)])
    if max_confs is not None:
        for command in commands:
            command.extend(["--max_confs", str(max_confs)])
    if patience is not None:
        for command in commands:
            command.extend(["--patience", str(patience)])
    if threshold is not None:
        for command in commands:
            command.extend(["--threshold", str(threshold)])
    if opt_steps is not None:
        for command in commands:
            command.extend(["--opt_steps", str(opt_steps)])
    for command in commands:
        command.extend(["--use_gpu", "True" if use_gpu else "False"])
    return commands


def _auto3d_env(mpi_np: int | None, cpu_workers: int | None) -> dict[str, str] | None:
    if mpi_np is None:
        return None
    per_process_threads = mpi_np
    if cpu_workers is not None and cpu_workers > 1:
        per_process_threads = max(1, mpi_np // cpu_workers)
    value = str(per_process_threads)
    return {
        "OMP_NUM_THREADS": value,
        "MKL_NUM_THREADS": value,
        "OPENBLAS_NUM_THREADS": value,
        "NUMEXPR_NUM_THREADS": value,
    }


def _cpu_worker_indices(cpu_workers: int) -> str:
    return ",".join(str(index) for index in range(cpu_workers))


def _find_output_sdf(output_dir: Path) -> Path | None:
    candidates = sorted(
        list(output_dir.glob("**/*_out.sdf")) + list(output_dir.glob("**/*_3d.sdf"))
    )
    return candidates[0] if candidates else None


def _python_executable() -> str:
    return shutil.which("python") or shutil.which("python3") or "python"


def _output_dir_name(output_sdf: Path) -> str:
    return output_sdf.parent.name


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

import multiprocessing as _mp


class _FakeManager:
    def Queue(self, maxsize: int = 0):
        return _mp.Queue(maxsize)


_mp.Manager = _FakeManager

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
