from __future__ import annotations

from pathlib import Path

from dsvr.reporting.progress import ProgressRecorder
from dsvr.workflow.status import run_status


def test_progress_files_and_status_survive_partial_run(tmp_path: Path) -> None:
    recorder = ProgressRecorder(tmp_path)

    recorder.record("Input validation", "completed", generated_count=1, accepted_count=1)
    recorder.record("Tautomer enumeration", "started", generated_count=4)

    assert (tmp_path / "progress.json").exists()
    assert (tmp_path / "progress.jsonl").exists()
    assert (tmp_path / "stage_summary.csv").exists()

    status = run_status(tmp_path)

    assert status["last_stage"] == "Tautomer enumeration"
    assert status["last_status"] == "started"
    assert status["stage_counts"]["Input validation"] == 1
