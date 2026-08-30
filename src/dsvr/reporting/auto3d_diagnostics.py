"""Consolidated Auto3D failure diagnostics (fix-auto3d-integration §5).

One :class:`Auto3DFailureBook` per run directory keeps:

- **root-cause records** — the first occurrence of an Auto3D execution
  failure is stored once (bounded excerpt, stage, engine, error class) in
  ``auto3d_root_causes.jsonl`` and surfaced once through the run's warning
  observers; affected candidates afterwards carry only a short reference
  (``auto3d_failed:<CLASS> (ref <id>)``) so warning stores stay readable;
- **failure memory** — invocation patterns that fail for infrastructure
  reasons (multiprocessing context crash, CUDA absence, timeout, missing
  Auto3D binary) or because Auto3D's own validation rejected the engine are
  not blindly re-attempted for the remaining molecules of the stage;
- **bounded notes** — helpers keeping per-candidate warning text short.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from dsvr.runners.auto3d_runner import (
    CLASS_ENGINE_INCOMPATIBLE,
    CLASS_EXECUTION,
    Auto3DUnavailableError,
    _emit_auto3d_notice_once,
    classify_auto3d_failure,
)

CLASS_UNAVAILABLE = "UNAVAILABLE"

#: Classes that mean a given invocation pattern cannot succeed for the rest
#: of the run and must not be re-attempted.
TERMINAL_MEMORY_CLASSES = frozenset(
    {"INFRA_MULTIPROCESSING", "CUDA_UNAVAILABLE", CLASS_UNAVAILABLE, CLASS_ENGINE_INCOMPATIBLE}
)

#: Infra classes whose failure has nothing to do with the chosen engine: one
#: such failure blocks the whole stage for every engine (a SemLock crash or
#: missing CUDA driver would recur identically).
STAGE_WIDE_MEMORY_CLASSES = frozenset(
    {"INFRA_MULTIPROCESSING", "CUDA_UNAVAILABLE", CLASS_UNAVAILABLE}
)

MAX_ROOT_EXCERPT_CHARS = 1024
MAX_CANDIDATE_NOTE_CHARS = 300

_VOLATILE_PATTERNS = (
    re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"),
    re.compile(r"/[^\s:]+"),
)


@dataclass(frozen=True)
class RootCause:
    root_cause_id: str
    stage: str
    error_class: str
    engine: str
    excerpt: str
    count: int = 1

    @property
    def short_note(self) -> str:
        return f"auto3d_failed:{self.error_class} (ref {self.root_cause_id})"


@dataclass
class Auto3DFailureBook:
    """Per-run root-cause store and terminal-failure memory."""

    run_dir: Path
    root_causes: dict[str, RootCause] = field(default_factory=dict)

    @property
    def _jsonl_path(self) -> Path:
        return self.run_dir / "auto3d_root_causes.jsonl"

    def terminal_reference(self, stage: str, engine: str) -> RootCause | None:
        """Return the terminal root cause blocking ``(stage, engine)``, if any.

        Infra failures (SemLock, CUDA absence, missing binary) block the
        whole stage regardless of engine; timeouts and engine-incompatibility
        rejections block only the same engine within the stage.
        """

        for cause in self.root_causes.values():
            if cause.stage != stage or cause.error_class not in TERMINAL_MEMORY_CLASSES:
                continue
            if cause.error_class in STAGE_WIDE_MEMORY_CLASSES or cause.engine == engine:
                return cause
        return None

    def record_failure(self, stage: str, engine: str, error: BaseException | str) -> RootCause:
        """Classify, deduplicate, and persist a failure; return its root cause."""

        if isinstance(error, Auto3DUnavailableError):
            error_class = CLASS_UNAVAILABLE
            text = str(error)
        elif isinstance(error, BaseException):
            text = str(error)
            error_class = CLASS_UNAVAILABLE if "not installed" in text else classify_auto3d_failure(text)
        else:
            text = error
            error_class = classify_auto3d_failure(text)
        # Aggregated downstream errors that already reference a recorded root
        # cause belong to it — bump its count instead of registering a
        # duplicate wrapper cause.
        match = re.search(r"\(ref ([0-9a-f]{12})\)", text)
        if match and match.group(1) in self.root_causes:
            existing = self.root_causes[match.group(1)]
            updated = RootCause(
                root_cause_id=existing.root_cause_id,
                stage=existing.stage,
                error_class=existing.error_class,
                engine=existing.engine,
                excerpt=existing.excerpt,
                count=existing.count + 1,
            )
            self.root_causes[existing.root_cause_id] = updated
            return updated
        excerpt = _normalize_excerpt(text)
        root_cause_id = hashlib.sha1(
            f"{stage}|{error_class}|{excerpt}".encode()
        ).hexdigest()[:12]
        existing = self.root_causes.get(root_cause_id)
        if existing is not None:
            updated = RootCause(
                root_cause_id=existing.root_cause_id,
                stage=existing.stage,
                error_class=existing.error_class,
                engine=existing.engine,
                excerpt=existing.excerpt,
                count=existing.count + 1,
            )
            self.root_causes[root_cause_id] = updated
            return updated
        cause = RootCause(
            root_cause_id=root_cause_id,
            stage=stage,
            error_class=error_class,
            engine=engine,
            excerpt=excerpt,
        )
        self.root_causes[root_cause_id] = cause
        self._write_root_cause(cause)
        return cause

    def _write_root_cause(self, cause: RootCause) -> None:
        self.run_dir.mkdir(parents=True, exist_ok=True)
        record = {
            "root_cause_id": cause.root_cause_id,
            "stage": cause.stage,
            "error_class": cause.error_class,
            "engine": cause.engine,
            "excerpt": cause.excerpt,
            "occurrences": cause.count,
        }
        with self._jsonl_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")
        # The full record is emitted once through the run's warning channel;
        # every later affected candidate only references root_cause_id.
        _emit_auto3d_notice_once(
            f"root-cause-{cause.root_cause_id}",
            (
                f"Auto3D failure root cause {cause.root_cause_id} "
                f"[{cause.error_class}] stage={cause.stage} engine={cause.engine}: "
                f"{cause.excerpt[:MAX_CANDIDATE_NOTE_CHARS]}"
            ),
        )

    def note_for_exception(self, stage: str, engine: str, error: BaseException | str) -> str:
        """Short per-candidate note referencing this exception's root cause.

        Errortext that already carries an ``auto3d_failed:... (ref ...)``
        note (e.g. an aggregated multi-engine failure) is referenced as-is —
        the true root cause was already recorded when first observed.
        """

        text = str(error)
        match = re.search(r"auto3d_failed:[A-Z_]+ \(ref [0-9a-f]{12}\)", text)
        if match:
            return match.group(0)
        return self.record_failure(stage, engine, error).short_note


_BOOKS: dict[Path, Auto3DFailureBook] = {}


def failure_book_for(run_dir: Path) -> Auto3DFailureBook:
    """Return the per-process :class:`Auto3DFailureBook` for a run directory."""

    key = Path(run_dir).resolve()
    book = _BOOKS.get(key)
    if book is None:
        book = Auto3DFailureBook(run_dir=key)
        _BOOKS[key] = book
    return book


def bounded(text: str, max_chars: int = MAX_CANDIDATE_NOTE_CHARS) -> str:
    """Bound a per-candidate warning string (task 5.3)."""

    if len(text) <= max_chars:
        return text
    return text[: max_chars - 1] + "…"


def _normalize_excerpt(text: str) -> str:
    head = text.strip().splitlines()[0] if text.strip() else ""
    for pattern in _VOLATILE_PATTERNS:
        head = pattern.sub("<x>", head)
    head = " ".join(head.split())
    excerpt = (head + "\n" + text.strip()) if len(text.strip()) > len(head) else head
    excerpt = " ".join(excerpt.split())
    return excerpt[:MAX_ROOT_EXCERPT_CHARS]


__all__ = [
    "CLASS_EXECUTION",
    "CLASS_UNAVAILABLE",
    "Auto3DFailureBook",
    "RootCause",
    "bounded",
    "failure_book_for",
]
