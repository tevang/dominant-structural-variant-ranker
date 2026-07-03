from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from dsvr.config import (
    RunConfig,
    load_config,
    merge_cli_overrides,
    write_resolved_config,
)


def test_default_config_resolves_ph_window() -> None:
    config = RunConfig()

    assert config.chemistry.ph == 7.0
    assert config.chemistry.ph_low == 7.0
    assert config.chemistry.ph_high == 7.0
    assert config.chemistry.solvent == "water"
    assert config.seeding.method == "etkdg"
    assert config.refinement.censo_enabled is False


def test_load_config_from_yaml(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    path = tmp_path / "config.yaml"
    path.write_text(
        "\n".join(
            [
                "run_name: test",
                f"input_path: {input_path}",
                "input_format: smi",
                "output_dir: run",
                "chemistry:",
                "  ph: 7.4",
                "  solvent: dmso",
                "enumeration:",
                "  max_protomers_per_molecule: 2",
                "seeding:",
                "  method: both",
                "refinement:",
                "  censo_enabled: true",
            ]
        ),
        encoding="utf-8",
    )

    config = load_config(path)

    assert config.run_name == "test"
    assert config.input_path == input_path
    assert config.chemistry.ph == 7.4
    assert config.chemistry.ph_low == 7.4
    assert config.chemistry.solvent == "dmso"
    assert config.enumeration.max_protomers_per_molecule == 2
    assert config.seeding.method == "both"
    assert config.refinement.censo_enabled is True


def test_load_resolved_config_anchors_relative_output_to_run_dir(tmp_path: Path) -> None:
    run_dir = tmp_path / "runs" / "dsvr"
    run_dir.mkdir(parents=True)
    path = run_dir / "resolved_config.yaml"
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    path.write_text(
        "\n".join(
            [
                "run_name: resolved",
                f"input_path: {input_path}",
                "input_format: smi",
                "output_dir: runs/dsvr",
            ]
        ),
        encoding="utf-8",
    )

    config = load_config(path)

    assert config.output_dir == run_dir


def test_cli_override_merge_updates_nested_fields(tmp_path: Path) -> None:
    config = merge_cli_overrides(
        RunConfig(),
        input_path=tmp_path / "input.sdf",
        output_dir=tmp_path / "out",
        ph=8.1,
        solvent="methanol",
        seeding_method="auto3d",
        censo_enabled=True,
    )

    assert config.input_path == tmp_path / "input.sdf"
    assert config.output_dir == tmp_path / "out"
    assert config.chemistry.ph == 8.1
    assert config.chemistry.ph_low == 8.1
    assert config.chemistry.ph_high == 8.1
    assert config.chemistry.solvent == "methanol"
    assert config.seeding.method == "auto3d"
    assert config.refinement.censo_enabled is True


def test_write_resolved_config(tmp_path: Path) -> None:
    config = RunConfig(output_dir=tmp_path / "run")

    path = write_resolved_config(config)

    assert path == tmp_path / "run" / "resolved_config.yaml"
    resolved = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert resolved["chemistry"]["ph_low"] == 7.0
    assert resolved["seeding"]["method"] == "etkdg"


def test_production_configs_load() -> None:
    for path in [
        Path("configs/production_balanced.yaml"),
        Path("configs/production_conservative.yaml"),
        Path("configs/exhaustive_debug.yaml"),
        Path("configs/crest_disk_safe.yaml"),
    ]:
        config = load_config(path)
        assert config.run_name


def test_rejects_partial_ph_window() -> None:
    with pytest.raises(ValidationError, match="ph_low and ph_high"):
        RunConfig(chemistry={"ph_low": 6.8})


def test_rejects_inverted_ph_window() -> None:
    with pytest.raises(ValidationError, match="ph_low must be <= ph_high"):
        RunConfig(chemistry={"ph_low": 8.0, "ph_high": 7.0})


def test_warns_for_unknown_solvent_but_allows_it() -> None:
    with pytest.warns(UserWarning, match="not in DSVR's conservative known-solvent list"):
        config = RunConfig(chemistry={"solvent": "custom-solvent"})

    assert config.chemistry.solvent == "custom-solvent"


def test_rejects_nonpositive_enumeration_caps() -> None:
    with pytest.raises(ValidationError, match="enumeration caps must be positive"):
        RunConfig(enumeration={"max_tautomers_per_protomer": 0})
