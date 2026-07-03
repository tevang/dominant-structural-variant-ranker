from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from dsvr.cli import app


def test_cli_commands_have_help() -> None:
    runner = CliRunner()
    commands = [
        "doctor",
        "validate-input",
        "enumerate-protomers",
        "enumerate-tautomers",
        "enumerate-stereo",
        "seed-etkdg",
        "seed-auto3d",
        "run-crest",
        "thermo",
        "rank",
        "run",
        "summarize",
    ]

    root = runner.invoke(app, ["--help"])
    assert root.exit_code == 0, root.output
    assert "--verbose" in root.output
    assert "--log-level" in root.output

    for command in commands:
        result = runner.invoke(app, [command, "--help"])
        assert result.exit_code == 0, result.output
        assert "Usage:" in result.output


def test_cli_summarize_missing_manifest_fails_actionably(tmp_path: Path) -> None:
    result = CliRunner().invoke(app, ["summarize", str(tmp_path)])

    assert result.exit_code != 0
    assert "manifest.json" in result.output
    assert "dsvr run" in result.output


def test_cli_run_fast_smoke_writes_expected_outputs(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCC(C)CC(C)(CCN)O 65986444\n", encoding="utf-8")
    config_path = tmp_path / "fast.yaml"
    config_path.write_text(
        "\n".join(
            [
                "run_name: cli-fast-smoke",
                f"input_path: {input_path}",
                "input_format: auto",
                f"output_dir: {tmp_path / 'configured-out'}",
                "enumeration:",
                "  max_protomers_per_molecule: 1",
                "  max_tautomers_per_protomer: 1",
                "  max_stereoisomers_per_tautomer: 1",
                "seeding:",
                "  method: etkdg",
                "  rdkit_num_conformers: 1",
                "crest:",
                "  enabled: false",
                "thermo:",
                "  enabled: false",
                "  xtb_hessian: false",
                "  xtb_thermo: false",
            ]
        ),
        encoding="utf-8",
    )
    outdir = tmp_path / "run"

    result = CliRunner().invoke(
        app,
        [
            "run",
            str(input_path),
            "--config",
            str(config_path),
            "--out",
            str(outdir),
            "--overwrite",
            "--no-resume",
        ],
    )

    assert result.exit_code == 0, result.output
    for name in [
        "manifest.json",
        "resolved_config.yaml",
        "logs/workflow.log",
        "inputs.csv",
        "ranked_variants.csv",
        "ranked_variants.json",
        "ranked_variants.sdf",
        "report.md",
    ]:
        assert (outdir / name).exists(), name
