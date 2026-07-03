from pathlib import Path

import yaml
from typer.testing import CliRunner

from dsvr.cli import app


def test_cli_overrides_config_values(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    config_path = tmp_path / "config.yaml"
    output_dir = tmp_path / "out"
    config_path.write_text(
        "\n".join(
            [
                "run_name: cli-test",
                f"input_path: {input_path}",
                "input_format: auto",
                "output_dir: should-be-overridden",
                "chemistry:",
                "  ph: 7.0",
                "  solvent: water",
                "seeding:",
                "  method: etkdg",
                "refinement:",
                "  censo_enabled: false",
            ]
        ),
        encoding="utf-8",
    )

    result = CliRunner().invoke(
        app,
        [
            "run",
            str(input_path),
            "--config",
            str(config_path),
            "--outdir",
            str(output_dir),
            "--ph",
            "8.2",
            "--solvent",
            "dmso",
            "--seeding-method",
            "auto3d",
            "--enable-censo",
        ],
    )

    assert result.exit_code == 0, result.output
    resolved = yaml.safe_load((output_dir / "resolved_config.yaml").read_text(encoding="utf-8"))
    assert resolved["input_path"] == str(input_path)
    assert resolved["output_dir"] == str(output_dir)
    assert resolved["chemistry"]["ph"] == 8.2
    assert resolved["chemistry"]["ph_low"] == 8.2
    assert resolved["chemistry"]["solvent"] == "dmso"
    assert resolved["seeding"]["method"] == "auto3d"
    assert resolved["refinement"]["censo_enabled"] is True


def test_cli_run_uses_config_output_dir_when_out_not_supplied(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    configured_out = tmp_path / "configured-out"
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "\n".join(
            [
                "run_name: cli-config-output",
                f"input_path: {input_path}",
                "input_format: auto",
                f"output_dir: {configured_out}",
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

    runner = CliRunner()
    result = runner.invoke(
        app,
        ["run", str(input_path), "--config", str(config_path)],
    )

    assert result.exit_code == 0, result.output
    assert (configured_out / "resolved_config.yaml").exists()
    assert not (tmp_path / "runs" / "dsvr" / "resolved_config.yaml").exists()


def test_cli_validate_input_writes_report_and_invalid_csv(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\nnot_a_smiles bad\n", encoding="utf-8")
    report_path = tmp_path / "validation_report.json"

    result = CliRunner().invoke(
        app,
        [
            "validate-input",
            str(input_path),
            "--format",
            "auto",
            "--out",
            str(report_path),
        ],
    )

    assert result.exit_code == 0, result.output
    report = yaml.safe_load(report_path.read_text(encoding="utf-8"))
    assert report["valid_count"] == 1
    assert report["invalid_count"] == 1
    assert (tmp_path / "invalid_inputs.csv").exists()
