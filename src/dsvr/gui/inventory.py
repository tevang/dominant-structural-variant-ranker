"""Read the canonical run output inventory written by the workflow.

The workflow publishes ``run_outputs.csv`` (see ``_publish_top_level_run_outputs``
in ``dsvr.workflow.engine``) with one row per artifact, encoding its path, kind,
existence, size and record count. This is the single source of truth the GUI's
artifact browser and view routing are built on.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Artifact:
    """Metadata for a single run output artifact."""

    name: str
    path: str
    kind: str
    exists: bool
    size_bytes: int
    record_count: int | str
    target: str


def iter_artifacts(inventory_path: Path) -> list[Artifact]:
    """Parse ``run_outputs.csv`` into artifact records."""
    artifacts: list[Artifact] = []
    if not inventory_path.exists():
        return artifacts
    with inventory_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            artifacts.append(
                Artifact(
                    name=row.get("artifact", ""),
                    path=row.get("path", ""),
                    kind=row.get("kind", ""),
                    exists=_as_bool(row.get("exists")),
                    size_bytes=_as_int(row.get("size_bytes")),
                    record_count=_count(row.get("record_count")),
                    target=row.get("target", ""),
                )
            )
    return artifacts


class RunInventory:
    """A parsed inventory of a run directory."""

    def __init__(self, artifacts: list[Artifact]) -> None:
        self.artifacts = artifacts

    def by_kind(self, kind: str) -> list[Artifact]:
        return [a for a in self.artifacts if a.kind == kind]

    def get(self, name: str) -> Artifact | None:
        for artifact in self.artifacts:
            if artifact.name == name:
                return artifact
        return None

    def exists(self, name: str) -> bool:
        artifact = self.get(name)
        return artifact is not None and artifact.exists

    def present_names(self) -> list[str]:
        return [a.name for a in self.artifacts]

    def kinds(self) -> list[str]:
        return sorted({a.kind for a in self.artifacts})


def load_inventory(run_dir: Path | str) -> RunInventory:
    return RunInventory(iter_artifacts(Path(run_dir) / "run_outputs.csv"))


def _as_bool(value: str | None) -> bool:
    return (value or "").lower() in {"true", "1", "yes"}


def _as_int(value: str | None) -> int:
    try:
        return int(value or 0)
    except ValueError:
        return 0


def _count(value: str | None) -> int | str:
    if value is None or value == "":
        return ""
    try:
        return int(value)
    except ValueError:
        return value
