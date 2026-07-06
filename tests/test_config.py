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
    assert config.seeding.auto3d_allow_rdkit_fallback is True
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
        Path("configs/auto3d_entropy_protocol.yaml"),
        Path("configs/auto3d_entropy_smoke.yaml"),
        Path("configs/ligprep_like_default.yaml"),
        Path("configs/ligprep_like_conservative.yaml"),
        Path("configs/ligprep_like_aggressive.yaml"),
        Path("configs/physics_validation_optional.yaml"),
    ]:
        config = load_config(path)
        assert config.run_name


def test_auto3d_entropy_protocol_config_loads() -> None:
    config = load_config(Path("configs/auto3d_entropy_protocol.yaml"))

    assert config.protocol == "auto3d_entropy"
    assert config.seeding.method == "auto3d"
    assert config.seeding.auto3d_internal_tautomer_stereo_enum is True
    assert config.seeding.auto3d_mpi_np == 28
    assert config.seeding.auto3d_cpu_workers == 28
    assert config.seeding.auto3d_memory_gb == 1
    assert config.seeding.auto3d_capacity == 1
    assert config.seeding.auto3d_model == "AIMNET"
    assert config.seeding.auto3d_k == 1
    assert config.seeding.auto3d_max_confs == 1
    assert config.seeding.auto3d_patience == 40
    assert config.seeding.auto3d_threshold == 0.5
    assert config.seeding.auto3d_opt_steps == 400
    assert config.crest.enabled is False
    assert config.thermo.enabled is False
    assert config.thermo.population_scope == "all_approximate"


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


def test_ligprep_like_default_config_loads() -> None:
    config = load_config(Path("configs/ligprep_like_default.yaml"))

    assert config.workflow_mode == "ligprep_like"
    assert config.chemistry.ph == 7.0
    assert config.chemistry.ph_low == 7.0
    assert config.protonation.max_protomers_per_molecule == 4
    assert config.tautomer_filtering.tool == "auto3d"
    assert config.tautomer_filtering.tauto_engine == "rdkit"
    assert config.tautomer_filtering.optimizing_engine == "ANI2xt"
    assert config.tautomer_filtering.max_rdkit_transforms == 256
    assert config.tautomer_filtering.fallback_if_timeout == "keep_input_and_canonical"
    assert config.stereoisomer_filtering.max_stereoisomers_per_tautomer == 16
    assert config.final_3d.k == 1
    assert config.final_3d.one_conformer_per_variant is True
    assert config.optional_validation.crest_xtb_enabled is False
    assert config.optional_validation.censo_enabled is False
    assert config.optional_validation.xtb_thermo_enabled is False


def test_ligprep_like_variant_configs_load() -> None:
    conservative = load_config(Path("configs/ligprep_like_conservative.yaml"))
    aggressive = load_config(Path("configs/ligprep_like_aggressive.yaml"))

    assert conservative.workflow_mode == "ligprep_like"
    assert conservative.protonation.max_protomers_per_molecule == 3
    assert aggressive.workflow_mode == "ligprep_like"
    assert aggressive.protonation.max_protomers_per_molecule == 8


def test_physics_validation_optional_config_loads() -> None:
    config = load_config(Path("configs/physics_validation_optional.yaml"))

    assert config.workflow_mode == "physics_validation"
    assert config.optional_validation.crest_xtb_enabled is True
    assert config.optional_validation.xtb_thermo_enabled is True


def test_agent_disabled_by_default() -> None:
    config = RunConfig()

    assert config.agent.enabled is False
    assert config.agent.backend == "ollama_codex_cli"
    assert config.agent.command == "codex --oss -m qwen3.6:35b"
    assert "classify_failure" in config.agent.allowed_tasks


def test_rejects_invalid_ligprep_caps_and_timeouts() -> None:
    with pytest.raises(ValidationError, match="protonation limits"):
        RunConfig(protonation={"max_protomers_per_molecule": -1})
    with pytest.raises(ValidationError, match="tautomer filtering limits"):
        RunConfig(tautomer_filtering={"timeout_seconds_per_protomer": -1})
    with pytest.raises(ValidationError, match="stereoisomer filtering limits"):
        RunConfig(stereoisomer_filtering={"max_stereoisomers_per_tautomer": -1})
    with pytest.raises(ValidationError, match="final_3d limits"):
        RunConfig(final_3d={"k": 0})


def test_ligprep_requires_one_final_conformer_per_variant() -> None:
    with pytest.raises(ValidationError, match="one_conformer_per_variant"):
        RunConfig(workflow_mode="ligprep_like", final_3d={"one_conformer_per_variant": False})
