from __future__ import annotations

import csv
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path


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
    active_command: str | None = None
    message: str = ""
    elapsed_seconds: float = 0.0
    run_dir_size_mb: float = 0.0
    xyz_file_count: int = 0


class ProgressRecorder:
    def __init__(self, run_dir: Path) -> None:
        self.run_dir = run_dir
        self.started = time.monotonic()
        self.events: list[ProgressEvent] = []
        self.progress_json = run_dir / "progress.json"
        self.progress_jsonl = run_dir / "progress.jsonl"
        self.stage_summary_csv = run_dir / "stage_summary.csv"

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
            active_command=active_command,
            message=message,
            elapsed_seconds=time.monotonic() - self.started,
            run_dir_size_mb=_directory_size_mb(self.run_dir),
            xyz_file_count=_xyz_file_count(self.run_dir),
        )
        self.events.append(event)
        self._write(event)

    def _write(self, event: ProgressEvent) -> None:
        self.run_dir.mkdir(parents=True, exist_ok=True)
        payload = {
            "last_event": asdict(event),
            "stage_counts": _stage_counts(self.events),
            "events_recorded": len(self.events),
        }
        self.progress_json.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with self.progress_jsonl.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(asdict(event), sort_keys=True) + "\n")
        _write_stage_summary(self.stage_summary_csv, self.events)


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
                "rejected_count",
                "skipped_count",
                "elapsed_seconds",
                "run_dir_size_mb",
                "xyz_file_count",
                "message",
            ],
        )
        writer.writeheader()
        for event in latest.values():
            row = asdict(event)
            writer.writerow({key: row.get(key) for key in writer.fieldnames or []})


def _directory_size_mb(path: Path) -> float:
    if not path.exists():
        return 0.0
    total = sum(item.stat().st_size for item in path.rglob("*") if item.is_file())
    return total / (1024**2)


def _xyz_file_count(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for item in path.rglob("*.xyz") if item.is_file())
