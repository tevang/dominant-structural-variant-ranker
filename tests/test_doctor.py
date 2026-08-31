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
    from dsvr.runners import unipka_runner

    # default config: no container runtime on this host → unipka check present but unavailable
    monkeypatch.setattr(unipka_runner.shutil, "which", lambda _name: None)

    statuses = check_tools(output_dir=tmp_path / "out")

    names = {status.name for status in statuses}
    assert {"python", "rdkit", "unipka", "molscrub", "Auto3D", "xtb", "crest"}.issubset(names)
    required = {status.name for status in statuses if status.required}
    assert {"python", "rdkit", "unipka", "xtb", "crest", "output-directory"}.issubset(required)
    assert "molscrub" not in required
    assert "Auto3D" not in required
    unipka = next(status for status in statuses if status.name == "unipka")
    assert not unipka.available
    assert "container" in unipka.detail or "Apptainer" in unipka.detail or "Uni-Pka" in unipka.detail


def test_doctor_groups_tool_interfaces_under_one_usability_row(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(tool_check, "executable_version", lambda *args, **kwargs: "mock 1.0")
    monkeypatch.setattr(tool_check, "python_import_check", lambda name: (False, None))
    monkeypatch.setattr(
        tool_check,
        "which_executable",
        lambda name: "/usr/bin/scrub.py" if name == "scrub.py" else None,
    )

    statuses = check_tools(output_dir=tmp_path / "out")

    groups = [status for status in statuses if status.kind == "tool"]
    assert [group.name for group in groups] == ["unipka", "molscrub", "Auto3D", "psi4"]

    # The summary row reports whether the tool is usable via any interface.
    molscrub = next(status for status in groups if status.name == "molscrub")
    assert not molscrub.required  # optional alternative to the Uni-Pka default
    assert molscrub.available
    assert "CLI" in molscrub.detail

    auto3d = next(status for status in groups if status.name == "Auto3D")
    assert not auto3d.required
    assert not auto3d.available

    # Interface rows are informational alternatives, never individually required.
    interface_rows = [status for status in statuses if status.group == "molscrub"]
    assert [status.name for status in interface_rows] == [
        "molscrub (python module)",
        "molscrub (CLI)",
    ]
    assert not any(status.required for status in interface_rows)
    assert not interface_rows[0].available
    assert interface_rows[1].available

    # Interface rows immediately follow their summary row.
    names = [status.name for status in statuses]
    summary_index = names.index("molscrub")
    assert names[summary_index + 1 : summary_index + 3] == [
        "molscrub (python module)",
        "molscrub (CLI)",
    ]


def test_doctor_payload_required_group_satisfied_by_any_interface(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(tool_check, "executable_version", lambda *args, **kwargs: "mock 1.0")
    monkeypatch.setattr(tool_check, "python_import_check", lambda name: (False, None))
    monkeypatch.setattr(
        tool_check,
        "which_executable",
        lambda name: "/usr/bin/scrub.py" if name == "scrub.py" else None,
    )
    from dsvr.runners import unipka_runner

    # pretend docker exists with a pullable image name → Uni-Pka required check satisfied
    monkeypatch.setattr(unipka_runner.shutil, "which", lambda name: f"/bin/{name}")
    monkeypatch.setattr(
        unipka_runner,
        "check_unipka_image",
        lambda config: config.container,
    )

    payload = tool_check.doctor_payload(output_dir=tmp_path / "out")

    assert "molscrub" not in payload["required_missing"]
    assert "unipka" not in payload["required_missing"]
    assert {"xtb", "crest"}.issubset(set(payload["required_missing"]))


def test_doctor_payload_flags_missing_unipka_image(
    tmp_path: Path,
    monkeypatch,
) -> None:
    from dsvr.runners import unipka_runner

    monkeypatch.setattr(unipka_runner.shutil, "which", lambda name: f"/bin/{name}")
    monkeypatch.setattr(unipka_runner, "check_unipka_image", lambda config: None)

    payload = tool_check.doctor_payload(output_dir=tmp_path / "out")

    assert "unipka" in payload["required_missing"]
    unipka_row = next(c for c in payload["checks"] if c["name"] == "unipka")
    assert "Zenodo" in unipka_row["detail"]


def test_cli_doctor_json_writes_machine_readable_report(
    tmp_path: Path,
    monkeypatch,
) -> None:
    monkeypatch.setattr(cli, "check_tools", lambda output_dir, protonation=None: _mock_statuses())
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
        lambda output_dir, protonation=None: [
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
