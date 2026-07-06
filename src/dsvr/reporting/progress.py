from __future__ import annotations

import csv
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from rich.console import Console


@dataclass(frozen=True)
class ProgressEvent:
    stage: str
    status: str
    molecule_index: int | None = None
    molecule_total: int | None = None
    molecule_name: str | None = None
    generated_count: int = 0
    accepted_count: int = 0
    rejected_count: int = 0
    skipped_count: int = 0
    timeout_count: int = 0
    active_command: str | None = None
    message: str = ""
    elapsed_seconds: float = 0.0
    run_dir_size_mb: float = 0.0
    xyz_file_count: int = 0


@dataclass(frozen=True)
class DiagnosticEvent:
    level: str
    stage: str
    message: str
    elapsed_seconds: float
    run_dir_size_mb: float


class ProgressRecorder:
    def __init__(
        self,
        run_dir: Path,
        *,
        terminal: bool = False,
        planned_stages: list[str] | None = None,
    ) -> None:
        self.run_dir = run_dir
        self.started = time.monotonic()
        self.events: list[ProgressEvent] = []
        self.progress_json = run_dir / "progress.json"
        self.progress_jsonl = run_dir / "progress.jsonl"
        self.stage_summary_csv = run_dir / "stage_summary.csv"
        self.variant_counts_csv = run_dir / "variant_counts.csv"
        self.warnings_jsonl = run_dir / "warnings.jsonl"
        self.failures_jsonl = run_dir / "failures.jsonl"
        self.terminal = terminal
        self.console = Console(stderr=True)
        self.planned_stages = planned_stages or []
        self._header_printed = False
        self._announced_stages: set[str] = set()
        self._completed_stages: set[str] = set()
        self.warning_count = 0
        self.failure_count = 0

    def record(
        self,
        stage: str,
        status: str,
        *,
        molecule_index: int | None = None,
        molecule_total: int | None = None,
        molecule_name: str | None = None,
        generated_count: int = 0,
        accepted_count: int = 0,
        rejected_count: int = 0,
        skipped_count: int = 0,
        timeout_count: int = 0,
        active_command: str | None = None,
        message: str = "",
    ) -> None:
        event = ProgressEvent(
            stage=stage,
            status=status,
            molecule_index=molecule_index,
            molecule_total=molecule_total,
            molecule_name=molecule_name,
            generated_count=generated_count,
            accepted_count=accepted_count,
            rejected_count=rejected_count,
            skipped_count=skipped_count,
            timeout_count=timeout_count,
            active_command=active_command,
            message=message,
            elapsed_seconds=time.monotonic() - self.started,
            run_dir_size_mb=_directory_size_mb(self.run_dir),
            xyz_file_count=_xyz_file_count(self.run_dir),
        )
        self.events.append(event)
        self._write(event)
        if self.terminal:
            self._print_terminal_event(event)

    def warning(self, stage: str, message: str) -> None:
        self.warning_count += 1
        self._write_diagnostic(self.warnings_jsonl, "warning", stage, message)

    def failure(self, stage: str, message: str) -> None:
        self.failure_count += 1
        self._write_diagnostic(self.failures_jsonl, "failure", stage, message)

    def _write(self, event: ProgressEvent) -> None:
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.warnings_jsonl.touch(exist_ok=True)
        self.failures_jsonl.touch(exist_ok=True)
        payload = {
            "last_event": asdict(event),
            "stage_counts": _stage_counts(self.events),
            "events_recorded": len(self.events),
            "warning_count": self.warning_count,
            "failure_count": self.failure_count,
        }
        self.progress_json.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with self.progress_jsonl.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(asdict(event), sort_keys=True) + "\n")
        _write_stage_summary(self.stage_summary_csv, self.events)
        _write_variant_counts(self.variant_counts_csv, self.events)

    def _write_diagnostic(self, path: Path, level: str, stage: str, message: str) -> None:
        self.run_dir.mkdir(parents=True, exist_ok=True)
        event = DiagnosticEvent(
            level=level,
            stage=stage,
            message=message,
            elapsed_seconds=time.monotonic() - self.started,
            run_dir_size_mb=_directory_size_mb(self.run_dir),
        )
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(asdict(event), sort_keys=True) + "\n")

    def _print_terminal_event(self, event: ProgressEvent) -> None:
        if not self._header_printed:
            self._print_header()
        if event.status == "started" and event.stage not in self._announced_stages:
            self._announced_stages.add(event.stage)
            command = f" command={event.active_command}" if event.active_command else ""
            self.console.print(
                f"[bold cyan]Stage {self._stage_position(event.stage)}: "
                f"{event.stage}[/bold cyan] [yellow]{event.status}[/yellow] "
                f"elapsed={_format_elapsed(event.elapsed_seconds)} "
                f"size={event.run_dir_size_mb:.2f}MB{command}"
            )
            return
        if event.status in {"running", "progress"}:
            self.console.print(self._progress_line(event))
            return
        if event.status in {"completed", "skipped"}:
            self._completed_stages.add(event.stage)
            self.console.print(self._completion_line(event))

    def _print_header(self) -> None:
        self._header_printed = True
        if self.planned_stages:
            self.console.print(
                f"[bold]Workflow progress[/bold]: {len(self.planned_stages)} stages planned; "
                f"run_dir=[bold]{self.run_dir}[/bold]"
            )
            self.console.print(
                "[dim]"
                + " -> ".join(
                    f"{index}. {stage}" for index, stage in enumerate(self.planned_stages, start=1)
                )
                + "[/dim]"
            )
        else:
            self.console.print(
                f"[bold]Workflow progress[/bold]: run_dir=[bold]{self.run_dir}[/bold]"
            )

    def _stage_position(self, stage: str) -> str:
        if stage in self.planned_stages:
            return f"{self.planned_stages.index(stage) + 1}/{len(self.planned_stages)}"
        return f"{len(self._announced_stages) + 1}/{len(self.planned_stages) or '?'}"

    def _progress_line(self, event: ProgressEvent) -> str:
        percent = _percent(event.molecule_index, event.molecule_total)
        count = _count_summary(event)
        molecule = f" item={event.molecule_name}" if event.molecule_name else ""
        command = f" command={event.active_command}" if event.active_command else ""
        size = f" size={event.run_dir_size_mb:.2f}MB"
        if percent is None:
            return (
                f"[cyan]{event.stage}[/cyan] [yellow]{event.status}[/yellow] "
                f"elapsed={_format_elapsed(event.elapsed_seconds)}{molecule}{count}{command}{size}"
            )
        return (
            f"[cyan]{event.stage}[/cyan] {_bar(percent)} {percent:5.1f}% "
            f"({event.molecule_index}/{event.molecule_total}) "
            f"elapsed={_format_elapsed(event.elapsed_seconds)}{molecule}{count}{command}{size}"
        )

    def _completion_line(self, event: ProgressEvent) -> str:
        overall = self._overall_percent()
        remaining = max(0.0, 100.0 - overall)
        count = _count_summary(event)
        return (
            f"[green]{event.stage} {event.status}[/green] "
            f"overall={_bar(overall)} {overall:5.1f}% "
            f"remaining={remaining:5.1f}% elapsed={_format_elapsed(event.elapsed_seconds)}"
            f"{count} size={event.run_dir_size_mb:.2f}MB"
        )

    def _overall_percent(self) -> float:
        if not self.planned_stages:
            return 0.0
        completed = len(
            [stage for stage in self.planned_stages if stage in self._completed_stages]
        )
        return min(100.0, 100.0 * completed / len(self.planned_stages))


def _stage_counts(events: list[ProgressEvent]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for event in events:
        counts[event.stage] = counts.get(event.stage, 0) + 1
    return counts


def _write_stage_summary(path: Path, events: list[ProgressEvent]) -> None:
    latest: dict[str, ProgressEvent] = {}
    for event in events:
        latest[event.stage] = event
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "stage",
                "status",
                "generated_count",
                "accepted_count",
                "selected_count",
                "rejected_count",
                "skipped_count",
                "timeout_count",
                "elapsed_seconds",
                "run_dir_size_mb",
                "xyz_file_count",
                "molecule_index",
                "molecule_total",
                "molecule_name",
                "active_command",
                "message",
            ],
        )
        writer.writeheader()
        for event in latest.values():
            row = asdict(event)
            row["selected_count"] = row["accepted_count"]
            writer.writerow({key: row.get(key) for key in writer.fieldnames or []})


def _write_variant_counts(path: Path, events: list[ProgressEvent]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "stage",
                "status",
                "generated_count",
                "selected_count",
                "rejected_count",
                "skipped_count",
                "timeout_count",
                "rejection_reason",
            ],
        )
        writer.writeheader()
        for event in events:
            if not any((event.generated_count, event.accepted_count, event.rejected_count, event.skipped_count, event.timeout_count)):
                continue
            writer.writerow(
                {
                    "stage": event.stage,
                    "status": event.status,
                    "generated_count": event.generated_count,
                    "selected_count": event.accepted_count,
                    "rejected_count": event.rejected_count,
                    "skipped_count": event.skipped_count,
                    "timeout_count": event.timeout_count,
                    "rejection_reason": event.message,
                }
            )


def _directory_size_mb(path: Path) -> float:
    if not path.exists():
        return 0.0
    total = 0
    for item in path.rglob("*"):
        try:
            if item.is_file():
                total += item.stat().st_size
        except OSError:
            continue
    return total / (1024**2)


def _xyz_file_count(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for item in path.rglob("*.xyz") if item.is_file())


def _percent(index: int | None, total: int | None) -> float | None:
    if index is None or total is None or total <= 0:
        return None
    return min(100.0, max(0.0, 100.0 * index / total))


def _bar(percent: float, *, width: int = 24) -> str:
    filled = int(round(width * max(0.0, min(100.0, percent)) / 100.0))
    return "|" + "#" * filled + "-" * (width - filled) + "|"


def _format_elapsed(seconds: float) -> str:
    total = int(seconds)
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{secs:02d}s"
    if minutes:
        return f"{minutes}m{secs:02d}s"
    return f"{secs}s"


def _count_summary(event: ProgressEvent) -> str:
    parts: list[str] = []
    if event.generated_count:
        parts.append(f"generated={event.generated_count}")
    if event.accepted_count:
        parts.append(f"selected={event.accepted_count}")
    if event.rejected_count:
        parts.append(f"rejected={event.rejected_count}")
    if event.skipped_count:
        parts.append(f"skipped={event.skipped_count}")
    if event.timeout_count:
        parts.append(f"timeouts={event.timeout_count}")
    return " " + " ".join(parts) if parts else ""


def read_jsonl_tail(path: Path, *, limit: int = 5) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines()[-limit:]:
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return rows
