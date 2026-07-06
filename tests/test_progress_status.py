from __future__ import annotations

import csv
from pathlib import Path

from dsvr.reporting.progress import ProgressRecorder
from dsvr.workflow.status import run_status


def test_progress_files_and_status_survive_partial_run(tmp_path: Path) -> None:
    (tmp_path / "resolved_config.yaml").write_text("run_name: test\n", encoding="utf-8")
    recorder = ProgressRecorder(tmp_path)

    recorder.record("Input validation", "completed", generated_count=1, accepted_count=1)
    recorder.record(
        "Tautomer filtering",
        "started",
        generated_count=4,
        accepted_count=2,
        rejected_count=2,
        timeout_count=1,
        active_command="Auto3D",
        message="outside energy window",
    )
    recorder.warning("Tautomer filtering", "mock warning")
    recorder.failure("Tautomer filtering", "mock failure")

    assert (tmp_path / "progress.json").exists()
    assert (tmp_path / "progress.jsonl").exists()
    assert (tmp_path / "stage_summary.csv").exists()
    assert (tmp_path / "variant_counts.csv").exists()
    assert (tmp_path / "warnings.jsonl").exists()
    assert (tmp_path / "failures.jsonl").exists()

    with (tmp_path / "variant_counts.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[-1]["selected_count"] == "2"
    assert rows[-1]["rejected_count"] == "2"
    assert rows[-1]["timeout_count"] == "1"
    assert rows[-1]["rejection_reason"] == "outside energy window"

    status = run_status(tmp_path)

    assert status["last_stage"] == "Tautomer filtering"
    assert status["last_status"] == "started"
    assert status["active_stage"] == "Tautomer filtering"
    assert status["last_completed_stage"] == "Input validation"
    assert status["stage_counts"]["Input validation"] == 1
    assert status["latest_warnings"][-1]["message"] == "mock warning"
    assert status["latest_failures"][-1]["message"] == "mock failure"
    assert status["resume_possible"] is True
