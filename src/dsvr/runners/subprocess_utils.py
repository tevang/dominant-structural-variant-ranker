from __future__ import annotations

import importlib.metadata
import importlib.util
import json
import os
import re
import shutil
import subprocess
import tempfile
import threading
import time
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

FAILURE_MARKERS = (
    "error",
    "failed",
    "segmentation fault",
    "nan",
    "not converged",
    "license",
    "command not found",
)


class ExternalToolError(RuntimeError):
    """Raised when an external command fails, times out, or cannot be launched."""

    def __init__(self, message: str, *, metadata: dict[str, object] | None = None) -> None:
        super().__init__(message)
        self.metadata = metadata or {}


@dataclass(frozen=True)
class LogTailSummary:
    path: Path
    line_count: int
    failure_markers: list[str] = field(default_factory=list)
    last_lines: list[str] = field(default_factory=list)


def run_command(
    command: list[str],
    cwd: Path | None = None,
    timeout_s: int | None = None,
    *,
    log_dir: Path | None = None,
    command_name: str | None = None,
    show_progress: bool = False,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run an external command while streaming stdout/stderr and writing metadata.

    The returned object intentionally remains ``subprocess.CompletedProcess`` compatible
    because individual chemistry runners inspect ``returncode`` and ``stdout``.
    """

    started = datetime.now(timezone.utc)  # noqa: UP017
    started_monotonic = time.monotonic()
    working_dir = Path(cwd) if cwd is not None else Path.cwd()
    run_log_dir = _prepare_log_dir(log_dir, command_name, started)
    stdout_path = run_log_dir / "stdout.log"
    stderr_path = run_log_dir / "stderr.log"
    combined_path = run_log_dir / "combined.log"
    metadata_path = run_log_dir / "command.json"
    metadata: dict[str, object] = {
        "command": command,
        "cwd": str(working_dir),
        "started_at": started.isoformat(),
        "ended_at": None,
        "duration_seconds": None,
        "returncode": None,
        "timeout_seconds": timeout_s,
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
        "combined_log": str(combined_path),
        "failure_markers": [],
    }
    _write_metadata(metadata_path, metadata)

    stdout_chunks: list[str] = []
    stderr_chunks: list[str] = []
    console = Console(stderr=True)

    try:
        process = subprocess.Popen(
            command,
            cwd=working_dir,
            env=os.environ | env if env is not None else None,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=1,
        )
    except OSError as exc:
        ended = datetime.now(timezone.utc)  # noqa: UP017
        metadata.update(
            {
                "ended_at": ended.isoformat(),
                "duration_seconds": time.monotonic() - started_monotonic,
                "launch_error": f"{type(exc).__name__}: {exc}",
                "returncode": None,
            }
        )
        _write_metadata(metadata_path, metadata)
        raise ExternalToolError(
            f"Failed to launch external command {' '.join(command)!r}: {exc}",
            metadata=metadata,
        ) from exc

    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_handle, combined_path.open("w", encoding="utf-8") as combined_handle:
        threads = [
            threading.Thread(
                target=_stream_pipe,
                args=(process.stdout, stdout_handle, combined_handle, stdout_chunks, "stdout"),
                daemon=True,
            ),
            threading.Thread(
                target=_stream_pipe,
                args=(process.stderr, stderr_handle, combined_handle, stderr_chunks, "stderr"),
                daemon=True,
            ),
        ]
        for thread in threads:
            thread.start()

        timed_out = False
        if show_progress:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                TimeElapsedColumn(),
                console=console,
                transient=True,
            ) as progress:
                task = progress.add_task(f"Running {' '.join(command[:2])}", total=None)
                timed_out = _wait_for_process(process, timeout_s, progress=progress, task=task)
        else:
            timed_out = _wait_for_process(process, timeout_s)

        if timed_out:
            process.kill()
        for thread in threads:
            thread.join(timeout=5)

    returncode = process.returncode
    ended = datetime.now(timezone.utc)  # noqa: UP017
    duration = time.monotonic() - started_monotonic
    stdout = "".join(stdout_chunks)
    stderr = "".join(stderr_chunks)
    summary = summarize_log(combined_path)
    metadata.update(
        {
            "ended_at": ended.isoformat(),
            "duration_seconds": duration,
            "returncode": returncode,
            "timed_out": timed_out,
            "failure_markers": summary.failure_markers,
        }
    )
    _write_metadata(metadata_path, metadata)

    completed = subprocess.CompletedProcess(
        args=command,
        returncode=-9 if timed_out and returncode is None else int(returncode or 0),
        stdout=stdout,
        stderr=stderr,
    )

    if timed_out:
        raise ExternalToolError(
            f"External command timed out after {timeout_s} s: {' '.join(command)}",
            metadata=metadata,
        )
    if check and completed.returncode != 0:
        raise ExternalToolError(
            f"External command failed with exit code {completed.returncode}: {' '.join(command)}",
            metadata=metadata,
        )
    return completed


def tail_log(path: Path, *, last_n: int = 20) -> list[str]:
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return lines[-last_n:]


def summarize_log(path: Path, *, last_n: int = 20) -> LogTailSummary:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines() if path.exists() else []
    markers = detect_failure_markers(lines)
    return LogTailSummary(
        path=path,
        line_count=len(lines),
        failure_markers=markers,
        last_lines=lines[-last_n:],
    )


def monitor_log(
    path: Path,
    *,
    interval_seconds: float = 10.0,
    stop_event: threading.Event | None = None,
    console: Console | None = None,
) -> None:
    """Tail a log at intervals and emit compact status summaries."""

    active_console = console or Console(stderr=True)
    seen_lines = 0
    while stop_event is None or not stop_event.is_set():
        if path.exists():
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
            new_lines = lines[seen_lines:]
            seen_lines = len(lines)
            markers = detect_failure_markers(new_lines)
            if markers:
                active_console.print(
                    f"[yellow]Detected possible failure marker(s): {', '.join(markers)}[/yellow]"
                )
            elif new_lines:
                active_console.print(summarize_progress(new_lines))
        time.sleep(interval_seconds)


def detect_failure_markers(lines: Iterable[str]) -> list[str]:
    text = "\n".join(lines).lower()
    return [marker for marker in FAILURE_MARKERS if marker in text]


def summarize_progress(lines: Iterable[str], *, max_chars: int = 240) -> str:
    compact = " ".join(line.strip() for line in lines if line.strip())
    if not compact:
        return "No new log output."
    if len(compact) > max_chars:
        compact = compact[-max_chars:]
    return f"Log progress: {compact}"


def which_executable(executable: str) -> str | None:
    return shutil.which(executable)


def executable_version(
    executable: str,
    *,
    args: list[str] | None = None,
    timeout_s: int = 10,
) -> str | None:
    path = which_executable(executable)
    if path is None:
        return None
    command = [path, *(args or ["--version"])]
    try:
        with tempfile.TemporaryDirectory(prefix="dsvr_version_probe_") as tmpdir:
            completed = run_command(
                command,
                timeout_s=timeout_s,
                check=False,
                log_dir=Path(tmpdir),
                command_name=f"{executable}_version",
            )
    except ExternalToolError:
        return None
    output = (completed.stdout or completed.stderr).strip()
    return output.splitlines()[0] if output else None


def python_import_check(module_name: str) -> tuple[bool, str | None]:
    if importlib.util.find_spec(module_name) is None:
        return False, None
    try:
        version = importlib.metadata.version(module_name)
    except importlib.metadata.PackageNotFoundError:
        version = None
    return True, version


def meets_minimum_version(version: str | None, minimum: str | None) -> bool | None:
    if minimum is None or version is None:
        return None
    version_tuple = _version_tuple(version)
    minimum_tuple = _version_tuple(minimum)
    return version_tuple >= minimum_tuple


def _wait_for_process(
    process: subprocess.Popen[str],
    timeout_s: int | None,
    *,
    progress: Progress | None = None,
    task: int | None = None,
) -> bool:
    deadline = None if timeout_s is None else time.monotonic() + timeout_s
    while process.poll() is None:
        if progress is not None and task is not None:
            progress.update(task)
        if deadline is not None and time.monotonic() > deadline:
            return True
        time.sleep(0.1)
    return False


def _stream_pipe(
    pipe: Iterable[str] | None,
    stream_handle,
    combined_handle,
    chunks: list[str],
    stream_name: str,
) -> None:
    if pipe is None:
        return
    for line in pipe:
        chunks.append(line)
        stream_handle.write(line)
        stream_handle.flush()
        combined_handle.write(f"[{stream_name}] {line}")
        combined_handle.flush()


def _prepare_log_dir(
    log_dir: Path | None,
    command_name: str | None,
    started: datetime,
) -> Path:
    base_dir = log_dir or Path.cwd() / "logs" / "subprocess"
    safe_name = _safe_name(command_name or "command")
    run_dir = base_dir / f"{started.strftime('%Y%m%dT%H%M%SZ')}_{safe_name}"
    suffix = 1
    candidate = run_dir
    while candidate.exists():
        suffix += 1
        candidate = run_dir.with_name(f"{run_dir.name}_{suffix}")
    candidate.mkdir(parents=True, exist_ok=False)
    return candidate


def _write_metadata(path: Path, metadata: dict[str, object]) -> None:
    path.write_text(json.dumps(metadata, indent=2, sort_keys=True, default=str) + "\n")


def _safe_name(value: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return sanitized[:80] or "command"


def _version_tuple(value: str) -> tuple[int, ...]:
    parts = re.findall(r"\d+", value)
    return tuple(int(part) for part in parts[:4])
