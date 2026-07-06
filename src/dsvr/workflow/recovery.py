from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from enum import StrEnum
from pathlib import Path
from typing import Any, Literal


CheckpointStatus = Literal["started", "completed", "failed", "skipped"]


class FailureKind(StrEnum):
    INPUT_ERROR = "INPUT_ERROR"
    PROTOMER_GENERATION_ERROR = "PROTOMER_GENERATION_ERROR"
    TAUTOMER_TIMEOUT = "TAUTOMER_TIMEOUT"
    TAUTOMER_ENUMERATION_ERROR = "TAUTOMER_ENUMERATION_ERROR"
    AUTO3D_FAILURE = "AUTO3D_FAILURE"
    STEREO_TIMEOUT = "STEREO_TIMEOUT"
    STEREO_ENUMERATION_ERROR = "STEREO_ENUMERATION_ERROR"
    EMBEDDING_FAILURE = "EMBEDDING_FAILURE"
    OPTIONAL_VALIDATION_FAILURE = "OPTIONAL_VALIDATION_FAILURE"
    DISK_LIMIT = "DISK_LIMIT"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class FailureRecord:
    failure_kind: FailureKind
    stage: str
    item_id: str | None
    item_name: str | None
    message: str
    action: str
    failed_at: str


def classify_failure(exc: BaseException, *, stage: str = "") -> FailureKind:
    text = f"{type(exc).__name__}: {exc}".lower()
    stage_key = stage.lower()
    if "disk" in text or "no space left" in text or "quota" in text:
        return FailureKind.DISK_LIMIT
    if "input" in stage_key or "smiles parse" in text or "invalid input" in text:
        return FailureKind.INPUT_ERROR
    if "protomer" in stage_key or "molscrub" in text or "protonation" in text:
        return FailureKind.PROTOMER_GENERATION_ERROR
    if "tautomer" in stage_key:
        return (
            FailureKind.TAUTOMER_TIMEOUT
            if _is_timeout(text)
            else FailureKind.TAUTOMER_ENUMERATION_ERROR
        )
    if "stereo" in stage_key or "stereoisomer" in stage_key:
        if _is_timeout(text):
            return FailureKind.STEREO_TIMEOUT
        if "embed" in text or "embedding" in text:
            return FailureKind.EMBEDDING_FAILURE
        return FailureKind.STEREO_ENUMERATION_ERROR
    if "auto3d" in stage_key or "auto3d" in text:
        return FailureKind.AUTO3D_FAILURE
    if "optional validation" in stage_key or "validation" in stage_key:
        return FailureKind.OPTIONAL_VALIDATION_FAILURE
    if "embed" in text or "embedding" in text:
        return FailureKind.EMBEDDING_FAILURE
    return FailureKind.UNKNOWN


def safe_action_for_failure(kind: FailureKind, *, stage: str = "") -> str:
    if kind == FailureKind.AUTO3D_FAILURE:
        return "retry_auto3d_or_skip_failed_variant"
    if kind == FailureKind.TAUTOMER_TIMEOUT:
        return "reduce_tautomer_cap_timeout"
    if kind == FailureKind.STEREO_TIMEOUT:
        return "reduce_stereo_cap_timeout"
    if kind == FailureKind.PROTOMER_GENERATION_ERROR:
        return "keep_fallback_parent_state"
    if kind == FailureKind.OPTIONAL_VALIDATION_FAILURE:
        return "skip_failed_variant"
    if kind == FailureKind.DISK_LIMIT:
        return "stop_or_skip_per_fail_fast"
    if "molecule" in stage.lower():
        return "skip_failed_molecule"
    return "skip_failed_variant"


class WorkflowRecoveryRecorder:
    def __init__(self, run_dir: Path) -> None:
        self.run_dir = run_dir
        self.root = run_dir / "checkpoints"
        self.stage_dir = self.root / "stages"
        self.molecule_dir = self.root / "molecules"
        self.failure_jsonl = self.root / "failures.jsonl"

    def stage(
        self,
        stage: str,
        status: CheckpointStatus,
        *,
        message: str = "",
        details: dict[str, Any] | None = None,
    ) -> None:
        payload = {
            "stage": stage,
            "status": status,
            "message": message,
            "updated_at": _now(),
            "details": details or {},
        }
        self.stage_dir.mkdir(parents=True, exist_ok=True)
        _write_json(self.stage_dir / f"{_safe_id(stage)}.json", payload)

    def molecule(
        self,
        *,
        item_id: str,
        item_name: str,
        stage: str,
        status: CheckpointStatus,
        failure_kind: FailureKind | None = None,
        action: str = "",
        message: str = "",
        details: dict[str, Any] | None = None,
    ) -> None:
        payload = {
            "item_id": item_id,
            "item_name": item_name,
            "stage": stage,
            "status": status,
            "failure_kind": failure_kind.value if failure_kind else None,
            "action": action,
            "message": message,
            "updated_at": _now(),
            "details": details or {},
        }
        self.molecule_dir.mkdir(parents=True, exist_ok=True)
        _write_json(self.molecule_dir / f"{_safe_id(item_id)}.json", payload)

    def failure(
        self,
        *,
        stage: str,
        item_id: str | None,
        item_name: str | None,
        exc: BaseException,
        action: str | None = None,
    ) -> FailureRecord:
        kind = classify_failure(exc, stage=stage)
        record = FailureRecord(
            failure_kind=kind,
            stage=stage,
            item_id=item_id,
            item_name=item_name,
            message=str(exc),
            action=action or safe_action_for_failure(kind, stage=stage),
            failed_at=_now(),
        )
        self.root.mkdir(parents=True, exist_ok=True)
        with self.failure_jsonl.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(asdict(record), sort_keys=True, default=str) + "\n")
        if item_id:
            self.molecule(
                item_id=item_id,
                item_name=item_name or item_id,
                stage=stage,
                status="failed",
                failure_kind=kind,
                action=record.action,
                message=record.message,
            )
        return record

    def molecule_state(self, item_id: str) -> dict[str, Any]:
        path = self.molecule_dir / f"{_safe_id(item_id)}.json"
        if not path.exists():
            return {}
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return {}


def should_skip_item_state(
    recorder: WorkflowRecoveryRecorder,
    item_id: str,
    *,
    resume: bool,
    stage: str | None = None,
) -> bool:
    if not resume:
        return False
    state = recorder.molecule_state(item_id)
    if stage is not None and state.get("stage") != stage:
        return False
    return state.get("status") in {"failed", "skipped"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()  # noqa: UP017


def _safe_id(value: str) -> str:
    safe = "".join(char if char.isalnum() or char in "._-" else "_" for char in value)
    return safe[:160] or "item"


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _is_timeout(text: str) -> bool:
    return "timeout" in text or "timed out" in text or "deadline" in text
