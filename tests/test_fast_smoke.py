from pathlib import Path

import pandas as pd

from dsvr.config import load_config
from dsvr.workflow.engine import run_workflow


def test_fast_smoke_runs_rdkit_enumeration_and_etkdg_without_external_tools(
    tmp_path: Path,
) -> None:
    config = load_config(Path("configs/fast_smoke.yaml")).model_copy(
        update={
            "input_path": Path("examples/test_molecules_minimal.smi"),
            "output_dir": tmp_path / "fast_smoke",
            "overwrite": True,
        }
    )

    result = run_workflow(config)

    assert result.molecule_count == 2
    assert (config.output_dir / "enumeration" / "tautomers" / "done.json").exists()
    assert (config.output_dir / "enumeration" / "stereoisomers" / "done.json").exists()
    assert (config.output_dir / "seeding" / "done.json").exists()
    assert (config.output_dir / "final_variants.sdf").exists()
    assert (config.output_dir / "final_variants.csv").exists()
    ranked = pd.read_csv(config.output_dir / "ranked_variants.csv")
    assert not ranked.empty
