from __future__ import annotations

import json
import subprocess
from pathlib import Path

import yaml
from typer.testing import CliRunner

from dsvr.agent.action_menu import parse_agent_output
from dsvr.agent.bug_package import build_bug_package
from dsvr.cli import app
from dsvr.config import RunConfig


def test_agent_cli_refuses_when_disabled(tmp_path: Path) -> None:
    run_dir = _write_run_dir(tmp_path, agent_enabled=False)

    result = CliRunner().invoke(app, ["agent", "diagnose", str(run_dir), "--latest"])

    assert result.exit_code != 0
    assert "Local diagnostic agent is disabled" in result.output


def test_agent_diagnose_uses_mocked_local_codex(tmp_path: Path, monkeypatch) -> None:
    run_dir = _write_run_dir(tmp_path, agent_enabled=True)
    _write_failure_fixture(run_dir)

    monkeypatch.setattr("dsvr.agent.local_qwen.shutil.which", lambda _cmd: "/usr/bin/codex")

    def fake_run(command, **kwargs):
        assert command == ["codex", "--oss", "-m", "qwen3.6:35b"]
        assert "Allowed actions:" in kwargs["input"]
        return subprocess.CompletedProcess(
            command,
            0,
            stdout=(
                "ACTION: retry_auto3d_cpu\n"
                "- Failure is in Auto3D.\n"
                "- CPU retry is in the allowed menu.\n"
                "- This does not change chemistry assumptions.\n"
            ),
            stderr="",
        )

    monkeypatch.setattr("dsvr.agent.local_qwen.subprocess.run", fake_run)

    result = CliRunner().invoke(app, ["agent", "diagnose", str(run_dir), "--latest"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["failure_kind"] == "AUTO3D_FAILURE"
    assert payload["deterministic_action"] == "retry_auto3d_cpu"
    assert payload["agent_action"] == "retry_auto3d_cpu"
    assert payload["agent_output_valid"] is True
    assert (run_dir / "bug_package" / "failure.json").exists()
    assert (run_dir / "bug_package" / "last_200_lines.log").exists()


def test_agent_unavailable_returns_human_review_without_failing(tmp_path: Path, monkeypatch) -> None:
    run_dir = _write_run_dir(tmp_path, agent_enabled=True)
    _write_failure_fixture(run_dir)
    monkeypatch.setattr("dsvr.agent.local_qwen.shutil.which", lambda _cmd: None)

    result = CliRunner().invoke(app, ["agent", "suggest-retry", str(run_dir), "--latest"])

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["agent_available"] is False
    assert payload["agent_action"] == "request_human_review"
    assert payload["deterministic_action"] == "retry_auto3d_cpu"


def test_invalid_agent_output_falls_back_to_human_review() -> None:
    decision = parse_agent_output("ACTION: delete_final_outputs\n- not allowed\n")

    assert decision.action == "request_human_review"
    assert decision.valid is False


def test_bug_package_contains_required_files(tmp_path: Path) -> None:
    run_dir = _write_run_dir(tmp_path, agent_enabled=True)
    config = RunConfig.model_validate(
        yaml.safe_load((run_dir / "resolved_config.yaml").read_text(encoding="utf-8"))
    )
    _write_failure_fixture(run_dir)

    package = build_bug_package(run_dir, config)

    assert package.path == run_dir / "bug_package"
    for filename in (
        "failure.json",
        "command.json",
        "last_200_lines.log",
        "config_fragment.yaml",
        "molecule.smi_or_sdf",
        "stage_summary.csv",
    ):
        assert (package.path / filename).exists()
    assert "AUTO3D_FAILURE" in (package.path / "failure.json").read_text(encoding="utf-8")


def _write_run_dir(tmp_path: Path, *, agent_enabled: bool) -> Path:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    config = RunConfig(
        input_path=input_path,
        output_dir=run_dir,
        agent={"enabled": agent_enabled},
    )
    (run_dir / "resolved_config.yaml").write_text(
        yaml.safe_dump(config.model_dump(mode="json"), sort_keys=False),
        encoding="utf-8",
    )
    return run_dir


def _write_failure_fixture(run_dir: Path) -> None:
    checkpoint_dir = run_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True)
    (checkpoint_dir / "failures.jsonl").write_text(
        json.dumps(
            {
                "failure_kind": "AUTO3D_FAILURE",
                "stage": "3D seeding",
                "item_id": "mol-1",
                "item_name": "ethanol",
                "message": "Auto3D GPU failed",
                "action": "retry_auto3d_or_skip_failed_variant",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    logs = run_dir / "logs" / "auto3d"
    logs.mkdir(parents=True)
    stdout = logs / "stdout.log"
    stderr = logs / "stderr.log"
    stdout.write_text("Auto3D start\n", encoding="utf-8")
    stderr.write_text("CUDA out of memory\n", encoding="utf-8")
    (logs / "combined.log").write_text("Auto3D start\nCUDA out of memory\n", encoding="utf-8")
    (logs / "command.json").write_text(
        json.dumps(
            {
                "command": ["auto3d"],
                "return_code": 1,
                "stdout_log": str(stdout),
                "stderr_log": str(stderr),
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (run_dir / "stage_summary.csv").write_text(
        "stage,status,generated_count\n3D seeding,failed,0\n",
        encoding="utf-8",
    )
