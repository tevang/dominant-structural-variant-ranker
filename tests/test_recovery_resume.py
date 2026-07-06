from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import yaml
from typer.testing import CliRunner

from dsvr.cli import app
from dsvr.config import RunConfig
from dsvr.workflow.engine import _run_step_list
from dsvr.workflow.recovery import FailureKind, classify_failure


def test_failure_classification_for_mocked_failures() -> None:
    assert (
        classify_failure(TimeoutError("tautomer timed out"), stage="Tautomer enumeration")
        == FailureKind.TAUTOMER_TIMEOUT
    )
    assert (
        classify_failure(RuntimeError("Auto3D failed"), stage="3D seeding")
        == FailureKind.AUTO3D_FAILURE
    )
    assert (
        classify_failure(RuntimeError("embedding failed"), stage="Stereoisomer enumeration")
        == FailureKind.EMBEDDING_FAILURE
    )
    assert (
        classify_failure(OSError("No space left on device"), stage="CREST/xTB")
        == FailureKind.DISK_LIMIT
    )


def test_step_failure_skips_failed_item_without_fail_fast(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    config = RunConfig(input_path=input_path, output_dir=tmp_path / "run")
    states = []
    items = [
        SimpleNamespace(id="ok", molname="ok"),
        SimpleNamespace(id="bad", molname="bad"),
        SimpleNamespace(id="ok2", molname="ok2"),
    ]

    def run_item(item):
        if item.id == "bad":
            raise RuntimeError("mock stereo enumeration failed")
        return [item.id]

    outputs = _run_step_list(
        "stereochemistry",
        items,
        config,
        states,
        run_item,
        progress_stage="Stereoisomer enumeration",
    )

    assert outputs == ["ok", "ok2"]
    assert (config.output_dir / "checkpoints" / "molecules" / "bad.json").exists()
    bad_state = json.loads(
        (config.output_dir / "checkpoints" / "molecules" / "bad.json").read_text(encoding="utf-8")
    )
    assert bad_state["status"] == "failed"
    assert bad_state["failure_kind"] == "STEREO_ENUMERATION_ERROR"
    failures = (config.output_dir / "checkpoints" / "failures.jsonl").read_text(encoding="utf-8")
    assert "mock stereo enumeration failed" in failures


def test_step_failure_raises_with_fail_fast(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    config = RunConfig(
        input_path=input_path,
        output_dir=tmp_path / "run",
        error_handling={"fail_fast": True},
    )
    item = SimpleNamespace(id="bad", molname="bad")

    def run_item(_item):
        raise RuntimeError("mock tautomer failure")

    try:
        _run_step_list(
            "tautomers",
            [item],
            config,
            [],
            run_item,
            progress_stage="Tautomer enumeration",
        )
    except RuntimeError as exc:
        assert "mock tautomer failure" in str(exc)
    else:
        raise AssertionError("fail_fast should re-raise the mocked failure")


def test_step_failure_raises_when_variant_skip_disabled(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    config = RunConfig(
        input_path=input_path,
        output_dir=tmp_path / "run",
        error_handling={"skip_failed_variant": False},
    )
    item = SimpleNamespace(id="bad", molname="bad")

    def run_item(_item):
        raise RuntimeError("mock stereo failure")

    try:
        _run_step_list(
            "stereochemistry",
            [item],
            config,
            [],
            run_item,
            progress_stage="Stereoisomer enumeration",
        )
    except RuntimeError as exc:
        assert "mock stereo failure" in str(exc)
    else:
        raise AssertionError("skip_failed_variant=false should re-raise variant failures")


def test_resume_command_loads_resolved_config_and_forces_resume(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "resolved_config.yaml").write_text(
        yaml.safe_dump(
            {
                "input_path": str(input_path),
                "output_dir": str(tmp_path / "wrong"),
                "resume": False,
                "overwrite": True,
            }
        ),
        encoding="utf-8",
    )
    seen = {}

    def fake_run_workflow(config):
        seen["config"] = config
        return SimpleNamespace(outdir=config.output_dir, molecule_count=1)

    monkeypatch.setattr("dsvr.cli.run_workflow", fake_run_workflow)
    result = CliRunner().invoke(app, ["resume", str(run_dir)])

    assert result.exit_code == 0, result.output
    assert seen["config"].output_dir == run_dir
    assert seen["config"].resume is True
    assert seen["config"].overwrite is False
