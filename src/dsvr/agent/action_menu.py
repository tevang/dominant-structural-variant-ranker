from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Any

ALLOWED_ACTIONS: tuple[str, ...] = (
    "retry_auto3d_cpu",
    "retry_auto3d_smaller_batch",
    "reduce_tautomer_cap_and_retry",
    "reduce_stereo_cap_and_retry",
    "keep_parent_tautomer_fallback",
    "skip_variant",
    "skip_molecule",
    "request_human_review",
)
REQUEST_HUMAN_REVIEW = "request_human_review"

_ACTION_RE = re.compile(r"^\s*(?:action\s*[:=]\s*)?([a-z0-9_]+)\s*$", re.IGNORECASE)


@dataclass(frozen=True)
class AgentDecision:
    action: str = REQUEST_HUMAN_REVIEW
    reasons: list[str] = field(default_factory=list)
    config_tweak: dict[str, Any] | None = None
    raw_output: str = ""
    valid: bool = False


def deterministic_action_for_failure(failure_kind: str | None) -> str:
    kind = (failure_kind or "").upper()
    if kind == "AUTO3D_FAILURE":
        return "retry_auto3d_cpu"
    if kind == "TAUTOMER_TIMEOUT":
        return "reduce_tautomer_cap_and_retry"
    if kind == "STEREO_TIMEOUT":
        return "reduce_stereo_cap_and_retry"
    if kind == "PROTOMER_GENERATION_ERROR":
        return "keep_parent_tautomer_fallback"
    if kind in {"OPTIONAL_VALIDATION_FAILURE", "EMBEDDING_FAILURE"}:
        return "skip_variant"
    if kind in {"INPUT_ERROR", "DISK_LIMIT", "UNKNOWN"}:
        return REQUEST_HUMAN_REVIEW
    return REQUEST_HUMAN_REVIEW


def parse_agent_output(text: str) -> AgentDecision:
    raw = text.strip()
    if not raw:
        return AgentDecision(raw_output=text)

    json_decision = _parse_json_decision(raw, text)
    if json_decision is not None:
        return json_decision

    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    action = _extract_action(lines)
    if action not in ALLOWED_ACTIONS:
        return AgentDecision(raw_output=text)

    reasons = [
        line.lstrip("-* ").strip()
        for line in lines
        if line.startswith(("-", "*")) and line.lstrip("-* ").strip()
    ][:3]
    return AgentDecision(
        action=action,
        reasons=reasons,
        raw_output=text,
        valid=True,
    )


def _parse_json_decision(raw: str, original: str) -> AgentDecision | None:
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError:
        return None
    if not isinstance(payload, dict):
        return AgentDecision(raw_output=original)
    action = payload.get("action")
    if action not in ALLOWED_ACTIONS:
        return AgentDecision(raw_output=original)
    reasons = payload.get("reasons")
    if not isinstance(reasons, list):
        reasons = []
    config_tweak = payload.get("config_tweak")
    if config_tweak is not None and not isinstance(config_tweak, dict):
        config_tweak = None
    return AgentDecision(
        action=str(action),
        reasons=[str(reason) for reason in reasons[:3]],
        config_tweak=config_tweak,
        raw_output=original,
        valid=True,
    )


def _extract_action(lines: list[str]) -> str:
    for line in lines[:4]:
        match = _ACTION_RE.match(line)
        if match:
            return match.group(1).lower()
    for action in ALLOWED_ACTIONS:
        if lines and action in lines[0].lower():
            return action
    return REQUEST_HUMAN_REVIEW
