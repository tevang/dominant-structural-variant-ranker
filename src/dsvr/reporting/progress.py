from __future__ import annotations

import csv
import json
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from rich.console import Console

from dsvr.workflow.recovery import WorkflowRecoveryRecorder


def _now_iso() -> str:
    return datetime.now().replace(microsecond=0).isoformat()


def _stage_key(stage_name: str | None) -> str | None:
    if stage_name is None:
        return None
    return (
        stage_name.strip()
        .lower()
        .replace("/", "_")
        .replace("-", "_")
        .replace(" ", "_")
    )


@dataclass
class WorkflowProgressState:
    stage_name: str | None = None
    stage: str | None = None
    total: int = 0
    running: int = 0
    succeeded: int = 0
    failed: int = 0
    skipped: int = 0
    completed: int = 0
    waiting: int = 0
    current_item: str | None = None
    current_item_name: str | None = None
    elapsed_seconds: float = 0.0
    rate_per_minute: float | None = None
    active_tool: str | None = None
    message: str | None = None
    last_update: str | None = None


class WorkflowProgress:
    """Lightweight workflow progress state, terminal display, and file writer."""

    def __init__(
        self,
        total_molecules: int,
        enabled: bool = True,
        *,
        run_dir: Path | None = None,
        interval_seconds: float = 2.0,
        console: Console | None = None,
        terminal: bool | None = None,
    ) -> None:
        self.total_molecules = total_molecules
        self.enabled = enabled
        self.run_dir = run_dir
        self.interval_seconds = max(0.0, interval_seconds)
        self.console = console or Console(stderr=True)
        self.terminal = enabled if terminal is None else terminal
        self.started = time.monotonic()
        self.stage_started = self.started
        self.state = WorkflowProgressState(total=total_molecules, waiting=total_molecules)
        self.progress_json = run_dir / "progress.json" if run_dir is not None else None
        self.progress_jsonl = run_dir / "progress.jsonl" if run_dir is not None else None
        self.stage_summary_csv = run_dir / "stage_summary.csv" if run_dir is not None else None
        self._events: list[dict[str, Any]] = []
        self._item_status: dict[str, str] = {}
        self._running_items: set[str] = set()
        self._last_terminal_update = 0.0
        self._live_terminal = bool(self.console.is_terminal and sys.stdout.isatty())

    def start_stage(self, stage_name: str, total_items: int | None = None) -> None:
        total = self.total_molecules if total_items is None else total_items
        self.stage_started = time.monotonic()
        self.state = WorkflowProgressState(
            stage_name=stage_name,
            stage=_stage_key(stage_name),
            total=total,
            waiting=total,
            last_update=_now_iso(),
        )
        self._item_status = {}
        self._running_items = set()
        self._emit("start_stage")

    def start_item(self, molecule_id: str, molecule_name: str | None = None) -> None:
        self._running_items.add(molecule_id)
        self.state.current_item = molecule_id
        self.state.current_item_name = molecule_name or molecule_id
        self._recount()
        self._emit("start_item", molecule_id=molecule_id, molecule_name=molecule_name)

    def finish_item(
        self,
        molecule_id: str,
        status: str = "success",
        message: str | None = None,
    ) -> None:
        normalized = "success" if status in {"success", "completed"} else status
        if normalized == "failed":
            self.fail_item(molecule_id, message=message)
            return
        if normalized == "skipped":
            self.skip_item(molecule_id, message=message)
            return
        self._running_items.discard(molecule_id)
        self._item_status[molecule_id] = "success"
        self.state.current_item = molecule_id
        self.state.current_item_name = self.state.current_item_name or molecule_id
        self.state.message = message
        self._recount()
        self._emit(
            "finish_item",
            molecule_id=molecule_id,
            molecule_name=self.state.current_item_name,
            status="success",
            message=message,
            force_terminal=True,
        )

    def fail_item(self, molecule_id: str, message: str | None = None) -> None:
        self._running_items.discard(molecule_id)
        self._item_status[molecule_id] = "failed"
        self.state.current_item = molecule_id
        self.state.current_item_name = self.state.current_item_name or molecule_id
        self.state.message = message
        self._recount()
        self._emit(
            "fail_item",
            molecule_id=molecule_id,
            molecule_name=self.state.current_item_name,
            status="failed",
            message=message,
            force_terminal=True,
        )

    def skip_item(self, molecule_id: str, message: str | None = None) -> None:
        self._running_items.discard(molecule_id)
        self._item_status[molecule_id] = "skipped"
        self.state.current_item = molecule_id
        self.state.current_item_name = self.state.current_item_name or molecule_id
        self.state.message = message
        self._recount()
        self._emit(
            "skip_item",
            molecule_id=molecule_id,
            molecule_name=self.state.current_item_name,
            status="skipped",
            message=message,
            force_terminal=True,
        )

    def update_message(self, message: str) -> None:
        self.state.message = message
        self._emit("update_message", message=message)

    def finish_stage(self) -> None:
        self._running_items.clear()
        self._recount()
        self._emit("finish_stage", force_terminal=True)

    def close(self) -> None:
        self._emit("close", force_terminal=True)

    def _recount(self) -> None:
        self.state.running = len(self._running_items)
        statuses = self._item_status.values()
        self.state.succeeded = sum(1 for status in statuses if status == "success")
        statuses = self._item_status.values()
        self.state.failed = sum(1 for status in statuses if status == "failed")
        statuses = self._item_status.values()
        self.state.skipped = sum(1 for status in statuses if status == "skipped")
        self.state.completed = self.state.succeeded + self.state.failed + self.state.skipped
        self.state.waiting = max(0, self.state.total - self.state.completed - self.state.running)
        self._refresh_timing()

    def _refresh_timing(self) -> None:
        self.state.elapsed_seconds = time.monotonic() - self.stage_started
        self.state.last_update = _now_iso()
        elapsed_minutes = self.state.elapsed_seconds / 60.0
        self.state.rate_per_minute = (
            self.state.completed / elapsed_minutes if elapsed_minutes > 0 else None
        )

    def _emit(
        self,
        event: str,
        *,
        molecule_id: str | None = None,
        molecule_name: str | None = None,
        status: str | None = None,
        message: str | None = None,
        force_terminal: bool = False,
    ) -> None:
        self._refresh_timing()
        record = {
            "timestamp": self.state.last_update,
            "stage": self.state.stage,
            "stage_name": self.state.stage_name,
            "event": event,
            "molecule_id": molecule_id,
            "molecule_name": molecule_name,
            "status": status,
            "completed": self.state.completed,
            "succeeded": self.state.succeeded,
            "failed": self.state.failed,
            "skipped": self.state.skipped,
            "running": self.state.running,
            "waiting": self.state.waiting,
            "total": self.state.total,
            "elapsed_seconds": self.state.elapsed_seconds,
            "rate_per_minute": self.state.rate_per_minute,
            "active_tool": self.state.active_tool,
            "current_item": self.state.current_item,
            "message": message if message is not None else self.state.message,
        }
        self._events.append(record)
        if not self.enabled:
            return
        self._write_files(record)
        self._print_terminal(record, force=force_terminal)

    def _write_files(self, record: dict[str, Any]) -> None:
        if self.run_dir is None:
            return
        self.run_dir.mkdir(parents=True, exist_ok=True)
        state_payload = asdict(self.state)
        state_payload["last_event"] = record
        state_payload["events_recorded"] = len(self._events)
        if self.progress_json is not None:
            self.progress_json.write_text(
                json.dumps(state_payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        if self.progress_jsonl is not None:
            with self.progress_jsonl.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(record, sort_keys=True) + "\n")
        self._write_stage_summary()

    def _write_stage_summary(self) -> None:
        if self.stage_summary_csv is None:
            return
        latest_by_stage: dict[str, dict[str, Any]] = {}
        for event in self._events:
            stage = str(event.get("stage_name") or event.get("stage") or "")
            latest_by_stage[stage] = event
        with self.stage_summary_csv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "stage",
                    "status",
                    "completed",
                    "succeeded",
                    "failed",
                    "skipped",
                    "waiting",
                    "total",
                    "elapsed_seconds",
                    "current_item",
                    "active_tool",
                    "last_update",
                    "message",
                ],
            )
            writer.writeheader()
            for stage, event in latest_by_stage.items():
                writer.writerow(
                    {
                        "stage": stage,
                        "status": event.get("event"),
                        "completed": event.get("completed", 0),
                        "succeeded": event.get("succeeded", 0),
                        "failed": event.get("failed", 0),
                        "skipped": event.get("skipped", 0),
                        "waiting": event.get("waiting", 0),
                        "total": event.get("total", 0),
                        "elapsed_seconds": f"{float(event.get('elapsed_seconds') or 0.0):.3f}",
                        "current_item": event.get("current_item") or "",
                        "active_tool": event.get("active_tool") or "",
                        "last_update": event.get("timestamp") or "",
                        "message": event.get("message") or "",
                    }
                )

    def _print_terminal(self, record: dict[str, Any], *, force: bool = False) -> None:
        if not self.terminal:
            return
        now = time.monotonic()
        event = str(record.get("event") or "")
        if not force and event not in {"start_stage", "finish_stage"}:
            if now - self._last_terminal_update < self.interval_seconds:
                return
        self._last_terminal_update = now
        if self._live_terminal:
            self.console.print(self._summary_block())
        elif event in {"finish_item", "fail_item", "skip_item", "finish_stage", "close"}:
            self.console.print(self._summary_line(), markup=False)

    def _summary_block(self) -> str:
        lines = [
            f"[bold cyan]Stage:[/bold cyan] {self.state.stage_name or ''}",
            (
                f"Molecules: {self.state.completed} / {self.state.total} completed | "
                f"Succeeded: {self.state.succeeded} | Failed: {self.state.failed} | "
                f"Skipped: {self.state.skipped}"
            ),
            f"Running: {self.state.running} | Waiting: {self.state.waiting}",
            f"Current molecule: {self.state.current_item_name or self.state.current_item or ''}",
            f"Elapsed: {_format_elapsed_colon(self.state.elapsed_seconds)}",
        ]
        if self.state.rate_per_minute is not None:
            lines.append(f"Rate: {self.state.rate_per_minute:.2f} molecules/min")
        if self.state.active_tool:
            lines.append(f"Tool: {self.state.active_tool}")
        if self.state.message:
            lines.append(f"Message: {self.state.message}")
        return "\n".join(lines)

    def _summary_line(self) -> str:
        return (
            f"[{self.state.stage or 'workflow'}] {self.state.completed}/{self.state.total} "
            f"completed | succeeded={self.state.succeeded} failed={self.state.failed} "
            f"skipped={self.state.skipped} waiting={self.state.waiting} | "
            f"current={self.state.current_item_name or self.state.current_item or ''}"
        )


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
        progress_interval: float = 2.0,
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
        self.recovery = WorkflowRecoveryRecorder(run_dir)
        self.console = Console(stderr=True)
        self.planned_stages = planned_stages or []
        self._header_printed = False
        self._announced_stages: set[str] = set()
        self._completed_stages: set[str] = set()
        self.progress_interval = max(0.0, progress_interval)
        self._last_terminal_update = 0.0
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
        if status in {"started", "completed", "failed", "skipped"}:
            self.recovery.stage(stage, status, message=message, details=asdict(event))
        self._write(event)
        if self.terminal and self._should_print_terminal_event(event):
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
        record = _workflow_event_from_progress_event(event)
        payload = {
            **_workflow_state_from_progress_event(event),
            "last_event": record,
            "legacy_last_event": asdict(event),
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
            handle.write(json.dumps(record, sort_keys=True) + "\n")
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

    def _should_print_terminal_event(self, event: ProgressEvent) -> bool:
        if event.status not in {"running", "progress"}:
            return True
        now = time.monotonic()
        if now - self._last_terminal_update < self.progress_interval:
            return False
        self._last_terminal_update = now
        return True

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


def _workflow_state_from_progress_event(event: ProgressEvent) -> dict[str, Any]:
    total = event.molecule_total or event.generated_count or 0
    completed = event.molecule_index or 0
    running = 0
    skipped = event.skipped_count
    failed = 0
    if event.status == "started":
        completed = 0
        running = 0
    elif event.status == "running":
        completed = event.molecule_index or 0
    elif event.status == "completed":
        completed = event.molecule_total or event.generated_count or completed
    elif event.status == "skipped":
        completed = total
        skipped = skipped or total
    elif event.status == "failed":
        completed = event.molecule_index or completed
        failed = 1
    succeeded = max(0, completed - failed - skipped)
    waiting = max(0, total - completed - running)
    elapsed = event.elapsed_seconds
    rate = (completed / (elapsed / 60.0)) if completed and elapsed > 0 else None
    return {
        "stage": _stage_key(event.stage),
        "stage_name": event.stage,
        "total": total,
        "running": running,
        "succeeded": succeeded,
        "failed": failed,
        "skipped": skipped,
        "completed": completed,
        "waiting": waiting,
        "current_item": event.molecule_name,
        "current_item_name": event.molecule_name,
        "elapsed_seconds": elapsed,
        "rate_per_minute": rate,
        "active_tool": event.active_command,
        "message": event.message or None,
        "last_update": _now_iso(),
    }


def _workflow_event_from_progress_event(event: ProgressEvent) -> dict[str, Any]:
    state = _workflow_state_from_progress_event(event)
    event_name = {
        "started": "start_stage",
        "running": "finish_item" if event.molecule_index is not None else "update_message",
        "progress": "update_message",
        "completed": "finish_stage",
        "failed": "fail_item",
        "skipped": "skip_item" if event.molecule_name else "finish_stage",
    }.get(event.status, event.status)
    status = None
    if event.status == "running":
        status = "success"
    elif event.status in {"failed", "skipped"}:
        status = event.status
    return {
        "timestamp": state["last_update"],
        "stage": state["stage"],
        "stage_name": state["stage_name"],
        "event": event_name,
        "molecule_id": event.molecule_name,
        "molecule_name": event.molecule_name,
        "status": status,
        "completed": state["completed"],
        "succeeded": state["succeeded"],
        "failed": state["failed"],
        "skipped": state["skipped"],
        "running": state["running"],
        "waiting": state["waiting"],
        "total": state["total"],
        "elapsed_seconds": state["elapsed_seconds"],
        "rate_per_minute": state["rate_per_minute"],
        "active_tool": state["active_tool"],
        "current_item": state["current_item"],
        "message": state["message"],
    }


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
                "completed",
                "succeeded",
                "failed",
                "skipped",
                "waiting",
                "total",
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
            row.update(_workflow_state_from_progress_event(event))
            row["stage"] = event.stage
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
            has_counts = any(
                (
                    event.generated_count,
                    event.accepted_count,
                    event.rejected_count,
                    event.skipped_count,
                    event.timeout_count,
                )
            )
            if not has_counts:
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


def _format_elapsed_colon(seconds: float) -> str:
    total = int(seconds)
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


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
