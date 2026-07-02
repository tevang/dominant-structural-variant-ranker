from pathlib import Path

from dsvr.config import DsvrConfig, load_config


def test_default_config() -> None:
    config = DsvrConfig()

    assert config.workflow.ph == 7.0
    assert config.workflow.solvent == "water"


def test_load_config(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text("workflow:\n  ph: 7.4\n  max_conformers: 3\n", encoding="utf-8")

    config = load_config(path)

    assert config.workflow.ph == 7.4
    assert config.workflow.max_conformers == 3

