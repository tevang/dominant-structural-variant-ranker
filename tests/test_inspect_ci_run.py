from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace


def _load_script_module():
    script_path = Path(__file__).resolve().parents[1] / "scripts" / "inspect_ci_run.py"
    spec = spec_from_file_location("inspect_ci_run", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_inspect_ci_run_parses_actions_url_and_builds_gh_command(monkeypatch) -> None:
    module = _load_script_module()
    captured = {}

    def fake_run(command, env):
        captured["command"] = command
        captured["env"] = env
        return SimpleNamespace(returncode=0)

    monkeypatch.setenv("GITHUB_TOKEN", "token-123")
    monkeypatch.delenv("GH_TOKEN", raising=False)
    monkeypatch.setattr(module.subprocess, "run", fake_run)

    exit_code = module.main(
        [
            "https://github.com/tevang/dominant-structural-variant-ranker/actions/runs/123456789",
        ]
    )

    assert exit_code == 0
    assert captured["command"] == [
        "gh",
        "run",
        "view",
        "123456789",
        "--repo",
        "tevang/dominant-structural-variant-ranker",
        "--log-failed",
    ]
    assert captured["env"]["GH_TOKEN"] == "token-123"


def test_inspect_ci_run_requires_repo_for_numeric_run_id(monkeypatch) -> None:
    module = _load_script_module()
    monkeypatch.setattr(
        module.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=0),
    )

    try:
        module.main(["123456789"])
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("expected SystemExit")
