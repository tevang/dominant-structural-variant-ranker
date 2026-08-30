from __future__ import annotations

import subprocess
import sys

import pytest

from dsvr.runners import auto3d_runner
from dsvr.runners.auto3d_runner import Auto3DExecutionError


def test_run_auto3d_stops_after_terminal_oscillation_failure(monkeypatch, tmp_path):
    calls: list[list[str]] = []

    monkeypatch.setattr(auto3d_runner, "_find_executable", lambda: "auto3d")
    monkeypatch.setattr(auto3d_runner.importlib.util, "find_spec", lambda name: object())

    def fake_run_command(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(
            args=command,
            returncode=1,
            stdout=(
                "Optimization finished at step 218: Total 3D structures: 1 "
                "Converged: 0 Dropped(Oscillating): 1 Active: 0\n"
                "OSError: File error: Invalid input file "
                "/tmp/job/auto3d_protomer_input_encoded_out.sdf\n"
            ),
            stderr="",
        )

    monkeypatch.setattr(auto3d_runner, "run_command", fake_run_command)

    with pytest.raises(Auto3DExecutionError):
        auto3d_runner.run_auto3d(
            tmp_path / "input.smi",
            tmp_path,
            k=1,
            model="AIMNET",
            internal_tautomer_stereo_enum=True,
        )

    assert len(calls) == 1


def test_run_auto3d_retries_nonterminal_failure(monkeypatch, tmp_path):
    calls: list[list[str]] = []

    monkeypatch.setattr(auto3d_runner, "_find_executable", lambda: "auto3d")
    monkeypatch.setattr(auto3d_runner.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(auto3d_runner, "_auto3d_major_version", lambda: 2)

    def fake_run_command(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(
            args=command,
            returncode=1,
            stdout="temporary command-line failure",
            stderr="",
        )

    monkeypatch.setattr(auto3d_runner, "run_command", fake_run_command)

    with pytest.raises(Auto3DExecutionError):
        auto3d_runner.run_auto3d(
            tmp_path / "input.smi",
            tmp_path,
            k=1,
            model="AIMNET",
            internal_tautomer_stereo_enum=True,
        )

    assert len(calls) == 3


def test_run_auto3d_v3_tries_single_cli_invocation(monkeypatch, tmp_path):
    calls: list[list[str]] = []

    monkeypatch.setattr(auto3d_runner, "_find_executable", lambda: "auto3d")
    monkeypatch.setattr(auto3d_runner.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(auto3d_runner, "_auto3d_major_version", lambda: 3)

    def fake_run_command(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(
            args=command,
            returncode=1,
            stdout="temporary v3 failure",
            stderr="",
        )

    monkeypatch.setattr(auto3d_runner, "run_command", fake_run_command)

    with pytest.raises(Auto3DExecutionError):
        auto3d_runner.run_auto3d(
            tmp_path / "input.smi",
            tmp_path,
            k=1,
            model="AIMNET",
            internal_tautomer_stereo_enum=True,
        )

    assert len(calls) == 1
    assert calls[0][0] == sys.executable
    assert calls[0][1].endswith("_auto3d_v3_wrapper.py")
    assert calls[0][2] == "run"


def _spawn_put(queue, value):
    queue.put(value)


def test_install_fake_manager_queue_survives_spawn_worker():
    """Regression test for the fork/spawn SemLock crash.

    Auto3D shares manager queues with spawn-context worker processes. A
    queue created in the default fork context raises "RuntimeError: A
    SemLock created in a fork context is being shared with a process in a
    spawn context" when the worker touches it; the fake manager must
    therefore build queues from the spawn context.
    """

    import multiprocessing as mp
    import multiprocessing.context as mp_context

    from dsvr.runners._auto3d_mp_shim import install_fake_manager

    original_manager = mp.Manager
    original_base_manager = mp_context.BaseContext.Manager
    try:
        install_fake_manager()
        queue = mp.Manager().Queue()
        process = mp.get_context("spawn").Process(
            target=_spawn_put, args=(queue, "sentinel")
        )
        process.start()
        process.join(timeout=60)
        assert process.exitcode == 0
        assert queue.get(timeout=30) == "sentinel"
    finally:
        mp.Manager = original_manager
        mp_context.BaseContext.Manager = original_base_manager


def test_run_auto3d_accepts_partial_output_on_nonzero_exit(monkeypatch, tmp_path):
    """Regression test for discarding partial Auto3D results.

    Auto3D exits nonzero when some inputs produce no output (e.g. charged
    molecules), but still writes an output SDF with the results it computed.
    run_auto3d must return that output so downstream stages can fill the
    missing variants, instead of discarding everything and burning time on
    duplicate GPU/CPU passes and per-molecule fallback retries.
    """

    monkeypatch.setattr(auto3d_runner, "_find_executable", lambda: "auto3d")
    monkeypatch.setattr(auto3d_runner.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(auto3d_runner, "_auto3d_major_version", lambda: 3)

    def fake_run_command(command, **kwargs):
        sdf = tmp_path / "job" / "input_out.sdf"
        sdf.parent.mkdir(parents=True, exist_ok=True)
        sdf.write_text("mol block\n$$$$\n", encoding="utf-8")
        return subprocess.CompletedProcess(
            args=command,
            returncode=6,
            stdout="",
            stderr="2 input molecule(s) produced no output",
        )

    monkeypatch.setattr(auto3d_runner, "run_command", fake_run_command)

    output, _command = auto3d_runner.run_auto3d(
        tmp_path / "input.smi",
        tmp_path,
        k=1,
        model="AIMNET",
        internal_tautomer_stereo_enum=False,
    )
    assert output.name == "input_out.sdf"


def test_run_auto3d_finds_output_next_to_input_file(monkeypatch, tmp_path):
    """Regression test: Auto3D creates job directories next to the INPUT file
    (<input_stem>_<job_name>/), which may be outside output_dir (tautomer
    filtering passes a sibling dir). run_auto3d must find that output."""

    input_dir = tmp_path / "input_here"
    output_dir = tmp_path / "output_elsewhere"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    input_smi = input_dir / "tautomer_candidates.smi"
    input_smi.write_text("CCO molA\n", encoding="utf-8")

    monkeypatch.setattr(auto3d_runner, "_find_executable", lambda: "auto3d")
    monkeypatch.setattr(auto3d_runner.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(auto3d_runner, "_auto3d_major_version", lambda: 3)

    def fake_run_command(command, **kwargs):
        job_name = command[command.index("--job-name") + 1]
        job_dir = input_dir / f"tautomer_candidates_{job_name}"
        job_dir.mkdir(parents=True, exist_ok=True)
        (job_dir / "tautomer_candidates_out.sdf").write_text(
            "mol block\n$$$$\n", encoding="utf-8"
        )
        return subprocess.CompletedProcess(args=command, returncode=0, stdout="", stderr="")

    monkeypatch.setattr(auto3d_runner, "run_command", fake_run_command)

    output, _command = auto3d_runner.run_auto3d(
        input_smi,
        output_dir,
        k=3,
        model="ANI2xt",
        internal_tautomer_stereo_enum=False,
    )
    assert output.name == "tautomer_candidates_out.sdf"


def test_find_output_sdf_prefers_current_job_over_stale_output(tmp_path):
    """Regression test for PR #3 review: a stale ``*_out.sdf`` left in a
    persistent output_dir by an earlier crashed invocation must not shadow
    this invocation's job-scoped output, even when it sorts first."""

    from dsvr.runners.auto3d_runner import _find_output_sdf

    stale = tmp_path / "aaa_stale" / "stale_out.sdf"
    stale.parent.mkdir(parents=True)
    stale.write_text("stale\n$$$$\n", encoding="utf-8")
    fresh = tmp_path / "zzz_input_abc123_x1" / "input_out.sdf"
    fresh.parent.mkdir(parents=True)
    fresh.write_text("mol block\n$$$$\n", encoding="utf-8")

    assert _find_output_sdf(tmp_path, job_name="abc123_x1") == fresh


def test_find_output_sdf_falls_back_to_unscoped_when_no_job_match(tmp_path):
    """Legacy Auto3D writes ``<input>_out.sdf`` without a job-named
    directory; the unscoped search must still find it."""

    from dsvr.runners.auto3d_runner import _find_output_sdf

    legacy = tmp_path / "input_out.sdf"
    legacy.write_text("mol block\n$$$$\n", encoding="utf-8")
    assert _find_output_sdf(tmp_path, job_name="nomatch") == legacy
