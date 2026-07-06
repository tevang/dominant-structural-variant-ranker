import json
from pathlib import Path

import pytest

from dsvr.runners import subprocess_utils
from dsvr.runners.subprocess_utils import (
    ExternalToolError,
    detect_failure_markers,
    run_command,
    summarize_log,
    summarize_progress,
)


def test_run_command_streams_logs_and_writes_metadata(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(subprocess_utils.subprocess, "Popen", _fake_popen(returncode=0))

    completed = run_command(
        ["mock-tool", "--flag"],
        cwd=tmp_path,
        log_dir=tmp_path / "logs",
        command_name="mock_tool",
    )

    assert completed.returncode == 0
    assert "normal progress" in completed.stdout
    run_dirs = list((tmp_path / "logs").iterdir())
    assert len(run_dirs) == 1
    assert (run_dirs[0] / "stdout.log").exists()
    assert (run_dirs[0] / "stderr.log").exists()
    assert (run_dirs[0] / "combined.log").exists()
    metadata = json.loads((run_dirs[0] / "command.json").read_text(encoding="utf-8"))
    assert metadata["command"] == ["mock-tool", "--flag"]
    assert metadata["returncode"] == 0
    assert metadata["return_code"] == 0
    assert metadata["elapsed_time_seconds"] >= 0
    assert metadata["stdout_log"].endswith("stdout.log")
    assert metadata["stderr_log"].endswith("stderr.log")
    assert metadata["ended_at"] is not None


def test_run_command_raises_external_tool_error_on_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        subprocess_utils.subprocess,
        "Popen",
        _fake_popen(returncode=2, stdout=["error: failed\n"]),
    )

    with pytest.raises(ExternalToolError) as excinfo:
        run_command(["mock-tool"], cwd=tmp_path, log_dir=tmp_path / "logs")

    assert "exit code 2" in str(excinfo.value)
    assert "error" in excinfo.value.metadata["failure_markers"]
    completed = run_command(
        ["mock-tool"],
        cwd=tmp_path,
        log_dir=tmp_path / "logs2",
        check=False,
    )
    assert completed.returncode == 2


def test_log_failure_marker_detection_and_progress_summary(tmp_path: Path) -> None:
    log_path = tmp_path / "combined.log"
    log_path.write_text("step 1\nnot converged\nsegmentation fault\n", encoding="utf-8")

    summary = summarize_log(log_path)

    assert "not converged" in summary.failure_markers
    assert "segmentation fault" in summary.failure_markers
    assert detect_failure_markers(["NaN in gradient"]) == ["nan"]
    assert summarize_progress(["iteration 1", "iteration 2"]).startswith("Log progress:")


def _fake_popen(
    *,
    returncode: int,
    stdout: list[str] | None = None,
    stderr: list[str] | None = None,
):
    class FakePopen:
        def __init__(self, *args, **kwargs) -> None:
            self.stdout = iter(stdout or ["normal progress\n"])
            self.stderr = iter(stderr or [])
            self.returncode = None
            self._final_returncode = returncode
            self._poll_count = 0

        def poll(self) -> int | None:
            self._poll_count += 1
            if self._poll_count > 1:
                self.returncode = self._final_returncode
                return self.returncode
            return None

        def kill(self) -> None:
            self.returncode = -9

    return FakePopen
