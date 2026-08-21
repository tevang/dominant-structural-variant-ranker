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
