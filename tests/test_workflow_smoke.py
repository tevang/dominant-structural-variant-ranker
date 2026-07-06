from pathlib import Path

from rdkit import Chem

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

    result = run_smoke_workflow(config=RunConfig(input_path=input_path, output_dir=outdir, protonation={"enabled": False}))

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


def test_workflow_auto3d_entropy_protocol_with_mocked_auto3d(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "auto3d-entropy"

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
        **kwargs,
    ) -> tuple[Path, list[str]]:
        assert internal_tautomer_stereo_enum is True
        # Auto3D protocol hands off a SMILES file: "SMILES ID" per line.
        first = None
        for line in input_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line:
                first = line
                break
        assert first is not None
        _, protomer_id = first.split(maxsplit=1)
        output_sdf = output_dir / "mock_auto3d_protocol.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for index, energy in enumerate((-10.0, -9.5), start=1):
            mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
            mol.SetProp("_Name", f"{protomer_id}_{index}")
            mol.SetProp("DSVR_PROTOMER_ID", protomer_id)
            mol.SetProp("E_kcal_mol", str(energy))
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "run", str(input_path), "--enumerate-tautomer"]

    def missing_molscrub(*args, **kwargs):
        raise engine_module.MolscrubUnavailableError("mock molscrub unavailable")

    monkeypatch.setattr(engine_module, "generate_protomer_candidates", missing_molscrub)
    monkeypatch.setattr(engine_module, "run_auto3d", fake_run_auto3d, raising=False)
    monkeypatch.setattr("dsvr.chemistry.conformers_auto3d.run_auto3d", fake_run_auto3d)

    result = run_workflow(
        RunConfig(
            protocol="auto3d_entropy",
            input_path=input_path,
            output_dir=outdir,
            protonation={"enabled": False},
            enumeration={"max_protomers_per_molecule": 1},
            seeding={"method": "auto3d", "auto3d_k": 2, "auto3d_model": "ANI2xt"},
            crest={"enabled": False},
            thermo={"enabled": False, "xtb_hessian": False, "xtb_thermo": False},
            variant_filtering={"enabled": False},
        )
    )

    assert result.molecule_count == 1
    assert (outdir / "seeding" / "auto3d_protocol" / "auto3d_protocol_seeds.sdf").exists()
    assert (outdir / "auto3d_representatives" / "auto3d_representative_scores.csv").exists()
    assert (outdir / "ranking" / "ranked_variants.csv").exists()
    assert (outdir / "ranked_variants.sdf").exists()
    assert (outdir / "auto3d_protocol_seeds.sdf").exists()
    assert (outdir / "auto3d_protocol_seeds.csv").exists()
    assert (outdir / "auto3d_adaptive_plan.csv").exists()
    assert (outdir / "all_protomers.sdf").exists()
    assert (outdir / "all_tautomers.sdf").exists()
    assert (outdir / "all_stereoisomers.sdf").exists()
    assert (outdir / "all_3d_conformers.sdf").exists()
    assert (outdir / "structure_generation_summary.csv").exists()
    assert (outdir / "structure_failures.csv").exists()
    assert (outdir / "run_outputs.csv").exists()
    summary = (outdir / "structure_generation_summary.csv").read_text(encoding="utf-8")
    assert "Auto3D representative generation" in summary
    assert "all_3d_conformers.sdf" in summary
    assert "ranked_variants.sdf" in summary
    manifest = (outdir / "manifest.json").read_text(encoding="utf-8")
    assert '"protocol": "auto3d_entropy"' in manifest


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
            protonation={"enabled": False},
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
        protonation={"enabled": False},
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
        protonation={"enabled": False},
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
