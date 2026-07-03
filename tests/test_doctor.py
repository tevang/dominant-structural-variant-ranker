import json
from pathlib import Path

from typer.testing import CliRunner

from dsvr import cli
from dsvr.models import ToolStatus
from dsvr.utils import tool_check
from dsvr.utils.tool_check import check_tools


def test_doctor_returns_default_workflow_and_optional_tool_statuses(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(tool_check, "executable_version", lambda *args, **kwargs: "mock 1.0")

    statuses = check_tools(output_dir=tmp_path / "out")

    names = {status.name for status in statuses}
    assert {"python", "rdkit", "molscrub", "Auto3D", "xtb", "crest"}.issubset(names)
    required = {status.name for status in statuses if status.required}
    assert {"python", "rdkit", "molscrub", "xtb", "crest", "output-directory"}.issubset(
        required
    )
    assert "Auto3D" not in required


def test_cli_doctor_json_writes_machine_readable_report(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(cli, "check_tools", lambda output_dir: _mock_statuses())
    json_out = tmp_path / "doctor.json"

    result = CliRunner().invoke(
        cli.app,
        [
            "doctor",
            "--output-dir",
            str(tmp_path / "out"),
            "--json",
            "--json-out",
            str(json_out),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(json_out.read_text(encoding="utf-8"))
    assert payload["ok"] is True
    assert payload["checks"][0]["name"] == "python"


def test_cli_doctor_strict_fails_only_for_required_missing(monkeypatch) -> None:
    monkeypatch.setattr(
        cli,
        "check_tools",
        lambda output_dir: [
            ToolStatus(
                name="xtb",
                kind="executable",
                required=True,
                available=False,
                detail="not on PATH",
            ),
            ToolStatus(
                name="Auto3D",
                kind="python-module",
                required=False,
                available=False,
                detail="optional",
            ),
        ],
    )

    non_strict = CliRunner().invoke(cli.app, ["doctor"])
    strict = CliRunner().invoke(cli.app, ["doctor", "--strict"])

    assert non_strict.exit_code == 0, non_strict.output
    assert strict.exit_code == 1
    assert "Required checks failed: xtb" in strict.output


def _mock_statuses() -> list[ToolStatus]:
    return [
        ToolStatus(
            name="python",
            kind="runtime",
            required=True,
            available=True,
            detail="mock-python",
            version="3.11.0",
        )
    ]
