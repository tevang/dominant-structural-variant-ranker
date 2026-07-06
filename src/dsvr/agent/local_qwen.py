from __future__ import annotations

import shlex
import shutil
import subprocess
from dataclasses import dataclass

from dsvr.agent.action_menu import (
    ALLOWED_ACTIONS,
    AgentDecision,
    parse_agent_output,
)
from dsvr.agent.policies import DISALLOWED_WITHOUT_EXPLICIT_USER_FLAG, sanitize_config_tweak
from dsvr.config import AgentConfig

SYSTEM_PROMPT = """You are a limited local diagnostic assistant.
Choose exactly one action from allowed action list.
Do not invent new actions.
Do not modify chemistry assumptions.
Explain in 3 bullet points why you chose action.
Return the first line as: ACTION: <allowed_action>.
"""


@dataclass(frozen=True)
class LocalAgentResult:
    available: bool
    decision: AgentDecision
    stdout: str = ""
    stderr: str = ""
    error: str = ""


def run_local_diagnostic_agent(
    *,
    agent: AgentConfig,
    task: str,
    bug_context: str,
) -> LocalAgentResult:
    args = shlex.split(agent.command)
    if not args or shutil.which(args[0]) is None:
        return LocalAgentResult(
            available=False,
            decision=AgentDecision(),
            error=f"Local agent command is unavailable: {agent.command}",
        )
    prompt = _build_prompt(task=task, bug_context=bug_context[: agent.max_context_chars])
    try:
        completed = subprocess.run(
            args,
            input=prompt,
            capture_output=True,
            text=True,
            timeout=agent.command_timeout_seconds,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return LocalAgentResult(
            available=False,
            decision=AgentDecision(),
            error=f"{type(exc).__name__}: {exc}",
        )
    if completed.returncode != 0:
        return LocalAgentResult(
            available=True,
            decision=AgentDecision(raw_output=completed.stdout),
            stdout=completed.stdout,
            stderr=completed.stderr,
            error=f"Local agent exited with code {completed.returncode}",
        )
    decision = parse_agent_output(completed.stdout)
    if decision.config_tweak is not None:
        decision = AgentDecision(
            action=decision.action,
            reasons=decision.reasons,
            config_tweak=sanitize_config_tweak(decision.config_tweak),
            raw_output=decision.raw_output,
            valid=decision.valid,
        )
    return LocalAgentResult(
        available=True,
        decision=decision,
        stdout=completed.stdout,
        stderr=completed.stderr,
    )


def _build_prompt(*, task: str, bug_context: str) -> str:
    return "\n".join(
        [
            SYSTEM_PROMPT.strip(),
            "",
            "Allowed actions:",
            "\n".join(f"- {action}" for action in ALLOWED_ACTIONS),
            "",
            "Disallowed actions without explicit user flag:",
            "\n".join(f"- {action}" for action in DISALLOWED_WITHOUT_EXPLICIT_USER_FLAG),
            "",
            f"Task: {task}",
            "",
            "Bug package:",
            bug_context,
            "",
        ]
    )
