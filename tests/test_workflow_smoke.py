from pathlib import Path

from dsvr.config import RunConfig
from dsvr.workflow.engine import run_smoke_workflow, run_workflow


def test_workflow_smoke_writes_outputs(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "run"

    result = run_smoke_workflow(config=RunConfig(input_path=input_path, output_dir=outdir))

    assert result.molecule_count == 1
    assert (outdir / "ranked.csv").exists()
    assert (outdir / "provenance.json").exists()
    assert (outdir / "resolved_config.yaml").exists()
    assert (outdir / "manifest.json").exists()
    assert (outdir / "logs" / "workflow.log").exists()
    assert (outdir / "ranking" / "ranked_variants.csv").exists()
    assert (outdir / "inputs.csv").exists()
    assert (outdir / "summary.md").exists()


def test_workflow_dry_run_writes_plan_without_external_execution(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "dry"

    result = run_workflow(
        RunConfig(input_path=input_path, output_dir=outdir, dry_run=True)
    )

    assert result.molecule_count == 0
    plan = outdir / "dry_run_plan.json"
    assert plan.exists()
    text = plan.read_text(encoding="utf-8")
    assert "protonation" in text
    assert "crest" in text
    assert (outdir / "manifest.json").exists()


def test_workflow_resume_skips_done_step_with_matching_hash(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "resume"
    config = RunConfig(
        input_path=input_path,
        output_dir=outdir,
        overwrite=False,
        resume=True,
        enumeration={
            "max_protomers_per_molecule": 1,
            "max_tautomers_per_protomer": 1,
            "max_stereoisomers_per_tautomer": 1,
        },
        seeding={"rdkit_num_conformers": 1},
        crest={"enabled": False},
        thermo={"enabled": False, "xtb_hessian": False, "xtb_thermo": False},
    )

    run_workflow(config)
    done = outdir / "enumeration" / "tautomers" / "done.json"
    first_done = done.read_text(encoding="utf-8")
    run_workflow(config)

    assert done.read_text(encoding="utf-8") == first_done
