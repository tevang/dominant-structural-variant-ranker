from pathlib import Path

from dsvr.config import DsvrConfig
from dsvr.workflow.engine import run_smoke_workflow


def test_workflow_smoke_writes_outputs(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "run"

    result = run_smoke_workflow(input_path=input_path, outdir=outdir, config=DsvrConfig())

    assert result.molecule_count == 1
    assert (outdir / "ranked.csv").exists()
    assert (outdir / "provenance.json").exists()
    assert (outdir / "summary.md").exists()

