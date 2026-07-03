from pathlib import Path

from dsvr.chemistry import tautomers as tautomer_module
from dsvr.chemistry.tautomers import TautomerEnumerationTimeout
from dsvr.config import RunConfig
from dsvr.workflow import engine as engine_module
from dsvr.workflow.engine import run_smoke_workflow, run_workflow
from dsvr.workflow.steps import mark_done, planned_steps, should_skip_step


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
    assert (outdir / "filtering" / "filtering_decisions.csv").exists()
    assert (outdir / "filtering" / "variant_penalties.csv").exists()
    assert (outdir / "filtering" / "accepted_variants.csv").exists()
    assert (outdir / "filtering" / "rejected_variants.csv").exists()
    assert (outdir / "filtering" / "penalty_breakdown.jsonl").exists()
    assert not (outdir / "seeding" / "rdkit" / "xyz").exists()


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


def test_workflow_tautomer_timeout_does_not_crash(tmp_path: Path, monkeypatch) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "timeout"

    def timeout_worker(*args, **kwargs):
        raise TautomerEnumerationTimeout("mock timeout")

    monkeypatch.setattr(tautomer_module, "_enumerate_molblocks_with_timeout", timeout_worker)

    result = run_workflow(
        RunConfig(
            input_path=input_path,
            output_dir=outdir,
            enumeration={"max_protomers_per_molecule": 1},
            seeding={"rdkit_num_conformers": 1},
            crest={"enabled": False},
            thermo={"enabled": False, "xtb_hessian": False, "xtb_thermo": False},
        )
    )

    assert result.molecule_count == 1
    report = (outdir / "report.md").read_text(encoding="utf-8")
    assert "Tautomer Timeout Counts" in report
    assert "tautomer enumeration timeout" in report


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


def test_workflow_resume_reuses_partial_tautomer_outputs(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "resume-partial"
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
    tautomer_done = outdir / "enumeration" / "tautomers" / "done.json"
    tautomer_done.unlink()

    def fail_if_recomputed(*args, **kwargs):
        raise AssertionError("existing tautomer output should be reused")

    monkeypatch.setattr(engine_module, "enumerate_tautomers", fail_if_recomputed)

    run_workflow(config)

    assert tautomer_done.exists()


def test_workflow_resume_rejects_failed_crest_checkpoint(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    config = RunConfig(input_path=input_path, output_dir=tmp_path / "run")
    step = {item.name: item for item in planned_steps(config)}["crest"]
    input_hash = "seed-record-hash"
    mark_done(step, input_hash, config, details={"count": 1, "enabled": True})
    provenance = step.output_dir / "mol" / "stereo" / "seed" / "crest_provenance.jsonl"
    provenance.parent.mkdir(parents=True)
    provenance.write_text(
        '{"crest_index": 0, "energy_kcal_mol": null}\n',
        encoding="utf-8",
    )

    assert not should_skip_step(step, input_hash, config)

    provenance.write_text(
        '{"crest_index": 1, "energy_kcal_mol": -42.0}\n',
        encoding="utf-8",
    )

    assert should_skip_step(step, input_hash, config)
