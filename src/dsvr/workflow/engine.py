from __future__ import annotations

import csv
import json
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.chemistry.conformers_auto3d import (
    generate_auto3d_seeds,
    generate_auto3d_seeds_from_protomers,
    score_auto3d_representative_variants,
)
from dsvr.chemistry.conformers_rdkit import generate_rdkit_seeds, read_stereo_sdf
from dsvr.chemistry.protonation import generate_protomer_candidates
from dsvr.chemistry.stereochemistry import enumerate_stereoisomers, read_tautomers_sdf
from dsvr.chemistry.tautomers import enumerate_tautomers, read_protomers_sdf
from dsvr.config import RunConfig, write_resolved_config
from dsvr.filtering.progress import decision_counts
from dsvr.filtering.selection import (
    FilteringDecision,
    select_seed_records,
    select_stereo_records,
    write_filtering_csv,
    write_filtering_decisions,
    write_penalty_outputs,
)
from dsvr.filtering.stereo_reduce import (
    StereoReductionResult,
    expand_enantiomer_mapped_crest_records,
    expand_enantiomer_mapped_thermo_records,
    reduce_seeds_for_crest,
    write_stereo_reduction_outputs,
)
from dsvr.filtering.xtb_prefilter import (
    XtbPrefilterDecision,
    apply_xtb_prefilter,
    write_xtb_prefilter_outputs,
)
from dsvr.io.read_inputs import read_molecules
from dsvr.io.write_outputs import write_final_ranked_outputs, write_json, write_ranked_csv
from dsvr.models import (
    AnyLineageRecord,
    CrestConformerRecord,
    MoleculeInput,
    ProtomerRecord,
    RankedVariantRecord,
    SeedConformerRecord,
    ThermoRecord,
    VariantRecord,
    WorkflowResult,
    make_protomer_id,
)
from dsvr.ranking.population import compute_delta_g_and_populations, write_ranked_outputs
from dsvr.reporting.audit import write_audit_tables
from dsvr.reporting.markdown import write_run_report, write_summary_markdown
from dsvr.reporting.progress import ProgressRecorder
from dsvr.runners.auto3d_runner import Auto3DExecutionError, Auto3DUnavailableError
from dsvr.runners.censo_runner import refine_top_ranked_with_censo
from dsvr.runners.crest_runner import read_seed_sdf, run_crest_for_seed
from dsvr.runners.molscrub_runner import MolscrubUnavailableError
from dsvr.runners.psi4_runner import rescore_top_ranked_with_psi4
from dsvr.runners.pyscf_runner import rescore_top_ranked_with_pyscf
from dsvr.runners.xtb_runner import run_xtb_thermo
from dsvr.utils.logging import configure_logging
from dsvr.workflow.provenance import build_provenance, write_all_provenance_outputs
from dsvr.workflow.steps import (
    StepState,
    file_hash,
    mark_done,
    planned_steps,
    records_hash,
    should_skip_step,
    skipped_state,
    write_dry_run_plan,
)


def run_workflow(config: RunConfig) -> WorkflowResult:
    outdir = config.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    configure_logging(level=config.logging.level, log_file=outdir / "logs" / "workflow.log")
    write_resolved_config(config, outdir)
    progress = ProgressRecorder(
        outdir,
        terminal=config.logging.level.upper() != "ERROR",
        planned_stages=_progress_stage_names(config),
    )
    progress.record("Input validation", "started")

    if config.dry_run:
        plan_path = write_dry_run_plan(config)
        manifest = _manifest(config, [], [], dry_run=True, extra={"dry_run_plan": str(plan_path)})
        write_json(outdir / "manifest.json", manifest)
        write_summary_markdown(outdir / "summary.md", molecule_count=0, variant_count=0)
        return WorkflowResult(outdir=outdir, molecule_count=0, ranked_records=[])

    steps = {step.name: step for step in planned_steps(config)}
    states: list[StepState] = []
    records: list[AnyLineageRecord] = []
    filtering_decisions: list[FilteringDecision] = []
    warnings: list[str] = []
    if config.variant_filtering.enabled and config.variant_filtering.mode == "exhaustive":
        warnings.append(
            "variant_filtering.mode=exhaustive is expensive; "
            "CREST/xTB may run on every generated seed."
        )
    input_hash = file_hash(config.input_path)

    molecules = read_molecules(
        config.input_path,
        input_format=config.input_format,
        invalid_output_path=outdir / "invalid_inputs.csv",
    )
    progress.record(
        "Input validation",
        "completed",
        generated_count=len(molecules),
        accepted_count=len(molecules),
    )
    _ensure_invalid_inputs_csv(outdir / "invalid_inputs.csv")
    _write_input_table(outdir / "input" / "inputs.csv", molecules)
    states.append(mark_done(steps["input"], input_hash, config, details={"count": len(molecules)}))

    if should_skip_step(steps["standardize"], records_hash(molecules), config):
        states.append(skipped_state(steps["standardize"], records_hash(molecules), config))
        molecules = _standardize_molecules(molecules, config)
        progress.record("Standardization", "skipped", generated_count=len(molecules))
    else:
        progress.record("Standardization", "started", generated_count=len(molecules))
        molecules = _standardize_molecules(molecules, config)
        _write_input_table(outdir / "input" / "standardized_inputs.csv", molecules)
        states.append(
            mark_done(
                steps["standardize"],
                records_hash(molecules),
                config,
                details={"enabled": config.chemistry.standardize},
            )
        )
        progress.record("Standardization", "completed", generated_count=len(molecules))

    protomers: list[ProtomerRecord] = []
    progress.record("Protomer generation", "started", generated_count=len(molecules))
    protomer_hash = records_hash(molecules)
    if should_skip_step(steps["protonation"], protomer_hash, config):
        states.append(skipped_state(steps["protonation"], protomer_hash, config))
        protomers = _load_protomers(outdir)
        progress.record("Protomer generation", "skipped", generated_count=len(protomers))
    else:
        for index, molecule in enumerate(molecules, start=1):
            existing = _load_existing_protomer_outputs(molecule, outdir, config)
            if existing:
                protomers.extend(existing)
            else:
                try:
                    protomers.extend(generate_protomer_candidates(molecule, config))
                except MolscrubUnavailableError as exc:
                    warnings.append(str(exc))
                    protomers.extend(_fallback_protomer(molecule, config))
            progress.record(
                "Protomer generation",
                "running",
                molecule_index=index,
                molecule_total=len(molecules),
                molecule_name=molecule.molname,
                generated_count=len(protomers),
            )
        states.append(
            mark_done(
                steps["protonation"],
                protomer_hash,
                config,
                details={"count": len(protomers), "fallback_warnings": warnings},
            )
        )
    records.extend(protomers)
    progress.record("Protomer generation", "completed", generated_count=len(protomers))

    if config.protocol == "auto3d_entropy":
        return _run_auto3d_entropy_protocol(
            config=config,
            molecules=molecules,
            protomers=protomers,
            records=records,
            states=states,
            warnings=warnings,
            progress=progress,
        )

    progress.record("Tautomer enumeration", "started", generated_count=len(protomers))
    tautomers = _run_step_list(
        "tautomers",
        protomers,
        config,
        states,
        lambda item: enumerate_tautomers(item, config),
        resume_loader=lambda: _load_tautomers(outdir),
        existing_loader=lambda item: _load_existing_tautomer_outputs(item, outdir, config),
        progress=progress,
        progress_stage="Tautomer enumeration",
    )
    records.extend(tautomers)
    progress.record("Tautomer enumeration", "completed", generated_count=len(tautomers))

    progress.record("Stereoisomer enumeration", "started", generated_count=len(tautomers))
    stereos = _run_step_list(
        "stereochemistry",
        tautomers,
        config,
        states,
        lambda item: enumerate_stereoisomers(item, config),
        resume_loader=lambda: _load_stereos(outdir),
        existing_loader=lambda item: _load_existing_stereo_outputs(item, outdir, config),
        progress=progress,
        progress_stage="Stereoisomer enumeration",
    )
    records.extend(stereos)
    progress.record("Stereoisomer enumeration", "completed", generated_count=len(stereos))
    progress.record("Cheap variant scoring", "started", generated_count=len(stereos))
    stereos_for_seeding, decisions = select_stereo_records(stereos, config, "pre_3d")
    filtering_decisions.extend(decisions)
    stereos_for_seeding, decisions = select_stereo_records(
        stereos_for_seeding,
        config,
        "cheap_score",
    )
    filtering_decisions.extend(decisions)
    progress.record(
        "Cheap variant scoring",
        "completed",
        generated_count=len(stereos),
        accepted_count=len(stereos_for_seeding),
        rejected_count=max(0, len(stereos) - len(stereos_for_seeding)),
    )

    seeds: list[SeedConformerRecord] = []
    progress.record("3D seeding", "started", generated_count=len(stereos_for_seeding))
    seed_input_hash = records_hash(stereos_for_seeding)
    if should_skip_step(steps["seeding"], seed_input_hash, config):
        states.append(skipped_state(steps["seeding"], seed_input_hash, config))
        seeds = _load_seeds(outdir)
        progress.record("3D seeding", "skipped", generated_count=len(seeds))
    else:
        seed_config = _config_for_seed_budget(config)
        if config.seeding.method in {"etkdg", "both"}:
            for index, stereo in enumerate(stereos_for_seeding, start=1):
                existing = _load_existing_seed_outputs(stereo, outdir, config)
                if existing:
                    seeds.extend(existing)
                else:
                    seeds.extend(generate_rdkit_seeds(stereo, seed_config))
                progress.record(
                    "3D seeding",
                    "running",
                    molecule_index=index,
                    molecule_total=len(stereos_for_seeding),
                    molecule_name=stereo.molname,
                    generated_count=len(seeds),
                )
        if config.seeding.method in {"auto3d", "both"}:
            try:
                progress.record(
                    "3D seeding",
                    "running",
                    generated_count=len(stereos_for_seeding),
                    message="Running Auto3D batch seeding.",
                )
                seeds.extend(generate_auto3d_seeds(stereos_for_seeding, seed_config))
            except (Auto3DExecutionError, Auto3DUnavailableError) as exc:
                warnings.append(str(exc))
                if config.seeding.method == "auto3d":
                    warnings.append(
                        "Auto3D failed in integrated workflow; falling back to RDKit ETKDG."
                    )
                    for index, stereo in enumerate(stereos_for_seeding, start=1):
                        seeds.extend(generate_rdkit_seeds(stereo, seed_config))
                        progress.record(
                            "3D seeding",
                            "running",
                            molecule_index=index,
                            molecule_total=len(stereos_for_seeding),
                            molecule_name=stereo.molname,
                            generated_count=len(seeds),
                        )
        states.append(
            mark_done(
                steps["seeding"],
                seed_input_hash,
                config,
                details={"count": len(seeds)},
            )
        )
    records.extend(seeds)
    progress.record("3D seeding", "completed", generated_count=len(seeds))
    progress.record("3D seed filtering", "started", generated_count=len(seeds))
    crest_seeds, decisions = select_seed_records(seeds, config)
    filtering_decisions.extend(decisions)
    progress.record("3D seed filtering", "completed", accepted_count=len(crest_seeds))
    progress.record("xTB prefilter", "started", generated_count=len(crest_seeds))
    crest_seeds, xtb_prefilter_decisions = apply_xtb_prefilter(crest_seeds, config)
    progress.record(
        "xTB prefilter",
        "completed",
        accepted_count=len(crest_seeds),
        rejected_count=len([item for item in xtb_prefilter_decisions if not item.selected]),
    )
    stereo_reduction = reduce_seeds_for_crest(crest_seeds, stereos_for_seeding, config)
    _write_filtering_outputs(
        outdir,
        filtering_decisions,
        stereo_reduction,
        xtb_prefilter_decisions,
    )

    progress.record("CREST/xTB", "started", generated_count=len(stereo_reduction.selected_seeds))
    representative_crest_records = _run_crest_or_seed_ranking(
        stereo_reduction.selected_seeds,
        config,
        states,
        warnings,
        progress,
    )
    progress.record(
        "CREST/xTB",
        "completed",
        generated_count=len(representative_crest_records),
        skipped_count=stereo_reduction.jobs_saved,
    )
    crest_records = expand_enantiomer_mapped_crest_records(
        representative_crest_records,
        {seed.id: seed for seed in crest_seeds},
        stereo_reduction,
        config,
    )
    records.extend(crest_records)

    thermo_records = []
    thermo_inputs = _select_thermo_inputs(representative_crest_records, config)
    progress.record("xTB thermo", "started", generated_count=len(thermo_inputs))
    thermo_input_hash = records_hash(thermo_inputs)
    if should_skip_step(steps["xtb_thermo"], thermo_input_hash, config):
        states.append(skipped_state(steps["xtb_thermo"], thermo_input_hash, config))
        thermo_records = _load_thermo_records(outdir)
        progress.record("xTB thermo", "skipped", generated_count=len(thermo_records))
    elif (
        config.thermo.enabled
        and thermo_inputs
        and _tool_available(config.crest.xtb_executable)
    ):
        for index, conformer in enumerate(thermo_inputs, start=1):
            thermo_records.append(run_xtb_thermo(conformer, config))
            progress.record(
                "xTB thermo",
                "running",
                molecule_index=index,
                molecule_total=len(thermo_inputs),
                molecule_name=conformer.molname,
                generated_count=len(thermo_records),
            )
        thermo_records = expand_enantiomer_mapped_thermo_records(
            thermo_records,
            crest_records,
            config,
        )
        states.append(
            mark_done(
                steps["xtb_thermo"],
                thermo_input_hash,
                config,
                details={
                    "count": len(thermo_records),
                    "input_count": len(thermo_inputs),
                    "enabled": config.thermo.enabled,
                    "max_variants_per_molecule": config.thermo.max_variants_per_molecule,
                    "max_conformers_per_variant": config.thermo.max_conformers_per_variant,
                },
            )
        )
    else:
        states.append(
            mark_done(
                steps["xtb_thermo"],
                thermo_input_hash,
                config,
                details={
                    "count": len(thermo_records),
                    "input_count": len(thermo_inputs),
                    "enabled": config.thermo.enabled,
                    "max_variants_per_molecule": config.thermo.max_variants_per_molecule,
                    "max_conformers_per_variant": config.thermo.max_conformers_per_variant,
                },
            )
        )
    records.extend(thermo_records)
    progress.record("xTB thermo", "completed", generated_count=len(thermo_records))

    rank_source = thermo_records if thermo_records else crest_records
    progress.record("Ranking", "started", generated_count=len(rank_source))
    ranked = compute_delta_g_and_populations(rank_source, config)
    write_ranked_outputs(ranked, outdir / "ranking")
    states.append(
        mark_done(
            steps["ranking"],
            records_hash(rank_source),
            config,
            details={"count": len(ranked)},
        )
    )
    records.extend(ranked)
    progress.record("Ranking", "completed", generated_count=len(ranked))

    censo_records = []
    progress.record("CENSO", "started", generated_count=len(ranked))
    if config.refinement.censo_enabled and _tool_available(config.refinement.censo_executable):
        censo_records = refine_top_ranked_with_censo(ranked, config)
    elif config.refinement.censo_enabled:
        warnings.append("CENSO requested but unavailable; integrated workflow skipped CENSO.")
    if censo_records:
        ranked_for_qm = censo_records
        records.extend(censo_records)
    else:
        ranked_for_qm = ranked
    states.append(
        mark_done(
            steps["censo"],
            records_hash(ranked),
            config,
            details={"count": len(censo_records), "enabled": config.refinement.censo_enabled},
        )
    )
    progress.record("CENSO", "completed", generated_count=len(censo_records))

    progress.record("QM rescoring", "started", generated_count=len(ranked_for_qm))
    qm_records = _run_optional_qm(ranked_for_qm, config)
    records.extend(qm_records)
    states.append(
        mark_done(
            steps["qm"],
            records_hash(ranked_for_qm),
            config,
            details={"count": len(qm_records), "backend": config.refinement.qm_backend},
        )
    )
    progress.record("QM rescoring", "completed", generated_count=len(qm_records))

    progress.record("Final reporting", "started", generated_count=len(ranked))
    write_all_provenance_outputs(records, outdir)
    write_audit_tables(
        outdir,
        records,
        filtering_decisions,
        xtb_prefilter_decisions,
        stereo_reduction,
    )
    write_final_ranked_outputs(outdir, ranked, config)
    write_summary_markdown(
        outdir / "summary.md",
        molecule_count=len(molecules),
        variant_count=len(ranked),
    )
    manifest = _manifest(config, states, warnings, dry_run=False)
    manifest["filtering"] = {
        "mode": config.variant_filtering.mode,
        "enabled": config.variant_filtering.enabled,
        "decision_counts": decision_counts(filtering_decisions),
        "decision_count": len(filtering_decisions),
        "stereo_reduction": {
            "jobs_saved": stereo_reduction.jobs_saved,
            "decision_count": len(stereo_reduction.decisions),
            "enabled": config.stereo_filtering.collapse_enantiomers_in_achiral_solvent
            and not config.stereo_filtering.solvent_is_chiral,
        },
        "xtb_prefilter": {
            "enabled": config.xtb_prefilter.enabled,
            "decision_count": len(xtb_prefilter_decisions),
            "pruned_count": len(
                [decision for decision in xtb_prefilter_decisions if not decision.selected]
            ),
        },
    }
    write_json(outdir / "manifest.json", manifest)
    write_run_report(
        outdir / "report.md",
        config=config,
        records=records,
        ranked_records=ranked,
        manifest=manifest,
        output_files=_final_output_files(outdir),
    )
    states.append(
        mark_done(
            steps["reports"],
            records_hash(records),
            config,
            details={"manifest": "manifest.json"},
        )
    )
    legacy_records = [
        VariantRecord(
            variant_id=record.id,
            parent_name=record.molname,
            smiles=record.isomeric_smiles,
            relative_energy_kcal_mol=record.relative_free_energy_kcal_mol,
            approximate_population=record.boltzmann_population,
            status="ranked",
        )
        for record in ranked
    ]
    write_ranked_csv(outdir / "ranked.csv", legacy_records)
    write_json(outdir / "provenance.json", build_provenance(config.input_path, config, molecules))
    progress.record("Final reporting", "completed", generated_count=len(ranked))
    return WorkflowResult(
        outdir=outdir,
        molecule_count=len(molecules),
        ranked_records=legacy_records,
    )


def run_smoke_workflow(config: RunConfig) -> WorkflowResult:
    return run_workflow(config)


def _progress_stage_names(config: RunConfig) -> list[str]:
    if config.protocol == "auto3d_entropy":
        return [
            "Input validation",
            "Standardization",
            "Protomer generation",
            "Auto3D representative generation",
            "Representative plausibility scoring",
            "Ranking",
            "Final reporting",
        ]
    return [
        "Input validation",
        "Standardization",
        "Protomer generation",
        "Tautomer enumeration",
        "Stereoisomer enumeration",
        "Cheap variant scoring",
        "3D seeding",
        "3D seed filtering",
        "xTB prefilter",
        "CREST/xTB",
        "xTB thermo",
        "Ranking",
        "CENSO",
        "QM rescoring",
        "Final reporting",
    ]


def _run_auto3d_entropy_protocol(
    *,
    config: RunConfig,
    molecules: list[MoleculeInput],
    protomers: list[ProtomerRecord],
    records: list[AnyLineageRecord],
    states: list[StepState],
    warnings: list[str],
    progress: ProgressRecorder,
) -> WorkflowResult:
    outdir = config.output_dir
    steps = {step.name: step for step in planned_steps(config)}

    for step_name in ("tautomers", "stereochemistry"):
        states.append(
            mark_done(
                steps[step_name],
                records_hash(protomers),
                config,
                details={
                    "count": 0,
                    "skipped_by_protocol": "auto3d_entropy",
                    "reason": "Auto3D performs tautomer/stereoisomer enumeration.",
                },
            )
        )

    progress.record("Auto3D representative generation", "started", generated_count=len(protomers))
    seed_hash = records_hash(protomers)
    if should_skip_step(steps["seeding"], seed_hash, config):
        states.append(skipped_state(steps["seeding"], seed_hash, config))
        seeds = _load_seeds(outdir)
        progress.record("Auto3D representative generation", "skipped", generated_count=len(seeds))
    else:
        progress.record(
            "Auto3D representative generation",
            "running",
            generated_count=len(protomers),
            message=(
                "Running Auto3D tautomer/stereo enumeration with "
                "representative conformer generation."
            ),
        )
        seeds = generate_auto3d_seeds_from_protomers(protomers, config)
        states.append(
            mark_done(
                steps["seeding"],
                seed_hash,
                config,
                details={
                    "count": len(seeds),
                    "protocol": "auto3d_entropy",
                    "auto3d_internal_tautomer_stereo_enum": True,
                },
            )
        )
    records.extend(seeds)
    progress.record("Auto3D representative generation", "completed", generated_count=len(seeds))

    empty_stereo_reduction = StereoReductionResult(
        selected_seeds=[],
        decisions=[],
        representative_to_equivalent_seed_ids={},
        jobs_saved=0,
    )
    _write_filtering_outputs(outdir, [], empty_stereo_reduction, [])

    for step_name in ("crest", "xtb_thermo"):
        states.append(
            mark_done(
                steps[step_name],
                records_hash(seeds),
                config,
                details={
                    "count": 0,
                    "enabled": False,
                    "skipped_by_protocol": "auto3d_entropy",
                },
            )
        )

    progress.record("Representative plausibility scoring", "started", generated_count=len(seeds))
    score_records = score_auto3d_representative_variants(seeds, config)
    records.extend(score_records)
    progress.record(
        "Representative plausibility scoring",
        "completed",
        generated_count=len(score_records),
    )

    progress.record("Ranking", "started", generated_count=len(score_records))
    ranked = compute_delta_g_and_populations(score_records, config)
    write_ranked_outputs(ranked, outdir / "ranking")
    states.append(
        mark_done(
            steps["ranking"],
            records_hash(score_records),
            config,
            details={"count": len(ranked), "source": "auto3d_representative_svp_score"},
        )
    )
    records.extend(ranked)
    progress.record("Ranking", "completed", generated_count=len(ranked))

    states.append(
        mark_done(
            steps["censo"],
            records_hash(ranked),
            config,
            details={"count": 0, "enabled": False, "skipped_by_protocol": "auto3d_entropy"},
        )
    )
    states.append(
        mark_done(
            steps["qm"],
            records_hash(ranked),
            config,
            details={"count": 0, "backend": "none", "skipped_by_protocol": "auto3d_entropy"},
        )
    )

    progress.record("Final reporting", "started", generated_count=len(ranked))
    write_all_provenance_outputs(records, outdir)
    write_audit_tables(outdir, records, [], [], empty_stereo_reduction)
    write_final_ranked_outputs(outdir, ranked, config)
    write_summary_markdown(
        outdir / "summary.md",
        molecule_count=len(molecules),
        variant_count=len(ranked),
    )
    manifest = _manifest(config, states, warnings, dry_run=False)
    manifest["protocol"] = "auto3d_entropy"
    manifest["filtering"] = {
        "mode": "protocol_disabled",
        "enabled": False,
        "decision_counts": {},
        "decision_count": 0,
        "stereo_reduction": {
            "jobs_saved": 0,
            "decision_count": 0,
            "enabled": False,
        },
        "xtb_prefilter": {
            "enabled": False,
            "decision_count": 0,
            "pruned_count": 0,
        },
    }
    write_json(outdir / "manifest.json", manifest)
    write_run_report(
        outdir / "report.md",
        config=config,
        records=records,
        ranked_records=ranked,
        manifest=manifest,
        output_files=_final_output_files(outdir),
    )
    states.append(
        mark_done(
            steps["reports"],
            records_hash(records),
            config,
            details={"manifest": "manifest.json"},
        )
    )
    legacy_records = [
        VariantRecord(
            variant_id=record.id,
            parent_name=record.molname,
            smiles=record.isomeric_smiles,
            relative_energy_kcal_mol=record.relative_free_energy_kcal_mol,
            approximate_population=record.boltzmann_population,
            status="ranked",
        )
        for record in ranked
    ]
    write_ranked_csv(outdir / "ranked.csv", legacy_records)
    write_json(outdir / "provenance.json", build_provenance(config.input_path, config, molecules))
    progress.record("Final reporting", "completed", generated_count=len(ranked))
    return WorkflowResult(
        outdir=outdir,
        molecule_count=len(molecules),
        ranked_records=legacy_records,
    )


def _run_step_list(
    step_name: str,
    inputs: list[Any],
    config: RunConfig,
    states: list[StepState],
    fn,
    resume_loader=None,
    existing_loader: Callable[[Any], list[Any]] | None = None,
    progress: ProgressRecorder | None = None,
    progress_stage: str | None = None,
) -> list[Any]:
    step = {item.name: item for item in planned_steps(config)}[step_name]
    input_hash = records_hash(inputs)
    if should_skip_step(step, input_hash, config):
        states.append(skipped_state(step, input_hash, config))
        if resume_loader is not None:
            loaded = resume_loader()
            if progress is not None and progress_stage is not None:
                progress.record(progress_stage, "skipped", generated_count=len(loaded))
            return loaded
        if progress is not None and progress_stage is not None:
            progress.record(progress_stage, "skipped")
        return []
    outputs = []
    for index, item in enumerate(inputs, start=1):
        existing = existing_loader(item) if existing_loader is not None else []
        if existing:
            outputs.extend(existing)
        else:
            outputs.extend(fn(item))
        if progress is not None and progress_stage is not None:
            progress.record(
                progress_stage,
                "running",
                molecule_index=index,
                molecule_total=len(inputs),
                molecule_name=getattr(item, "molname", None),
                generated_count=len(outputs),
            )
    states.append(mark_done(step, input_hash, config, details={"count": len(outputs)}))
    return outputs


def _config_for_seed_budget(config: RunConfig) -> RunConfig:
    if (
        not config.variant_filtering.enabled
        or config.variant_filtering.mode == "exhaustive"
    ):
        return config
    data = config.model_dump(mode="python")
    data["seeding"]["auto3d_k"] = min(
        config.seeding.auto3d_k,
        config.variant_filtering.max_seeds_per_variant,
    )
    return RunConfig.model_validate(data)


def _write_filtering_outputs(
    outdir: Path,
    decisions: list[FilteringDecision],
    stereo_reduction: StereoReductionResult | None = None,
    xtb_prefilter_decisions: list[XtbPrefilterDecision] | None = None,
) -> None:
    write_filtering_decisions(outdir / "filtering" / "filtering_decisions.jsonl", decisions)
    write_filtering_csv(outdir / "filtering" / "filtering_decisions.csv", decisions)
    write_penalty_outputs(outdir / "filtering", decisions)
    if stereo_reduction is not None:
        write_stereo_reduction_outputs(outdir / "filtering", stereo_reduction)
    if xtb_prefilter_decisions is not None:
        write_xtb_prefilter_outputs(outdir / "filtering", xtb_prefilter_decisions)


def _standardize_molecules(
    molecules: list[MoleculeInput],
    config: RunConfig,
) -> list[MoleculeInput]:
    if not config.chemistry.standardize:
        return molecules
    standardized = []
    for molecule in molecules:
        mol = Chem.Mol(molecule.rdkit_mol)
        Chem.SanitizeMol(mol)
        canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        standardized.append(
            molecule.model_copy(
                update={
                    "rdkit_mol": mol,
                    "canonical_smiles": canonical,
                    "isomeric_smiles": isomeric,
                }
            )
        )
    return standardized


def _fallback_protomer(molecule: MoleculeInput, config: RunConfig) -> list[ProtomerRecord]:
    output_dir = config.output_dir / "enumeration" / "protomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    mol = Chem.Mol(molecule.rdkit_mol)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    metadata = {
        "fallback": True,
        "reason": "molscrub unavailable; original molecule used as single protomer candidate",
    }
    record = ProtomerRecord(
        id=make_protomer_id(molecule.input_id, 1, canonical, isomeric, metadata),
        parent_id=molecule.input_id,
        input_molecule_id=molecule.input_id,
        molname=molecule.molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula=formula,
        formal_charge=Chem.GetFormalCharge(mol),
        explicit_proton_count=_explicit_proton_count(mol),
        source_software="dsvr-fallback",
        source_python_function="dsvr.workflow.engine._fallback_protomer",
        output_paths=[
            output_dir / f"{molecule.input_id}_protomers.sdf",
            output_dir / f"{molecule.input_id}_protomers.csv",
        ],
        warnings=[
            "molscrub unavailable; original molecule used as a single protomer candidate. "
            "This is suitable for smoke tests only."
        ],
        metadata=metadata,
        protomer_index=1,
        rdkit_mol=mol,
    )
    _write_protomer_outputs(record)
    return [record]


def _run_crest_or_seed_ranking(
    seeds: list[SeedConformerRecord],
    config: RunConfig,
    states: list[StepState],
    warnings: list[str],
    progress: ProgressRecorder | None = None,
) -> list[CrestConformerRecord]:
    step = {item.name: item for item in planned_steps(config)}["crest"]
    seed_hash = records_hash(seeds)
    if should_skip_step(step, seed_hash, config):
        states.append(skipped_state(step, seed_hash, config))
        records = _load_crest_records(config.output_dir)
        if progress is not None:
            progress.record("CREST/xTB", "skipped", generated_count=len(records))
        return records
    records = []
    crest_tools_available = _tool_available(config.crest.executable) and _tool_available(
        config.crest.xtb_executable
    )
    if config.crest.enabled and crest_tools_available:
        for index, seed in enumerate(seeds, start=1):
            records.extend(run_crest_for_seed(seed, config))
            if progress is not None:
                progress.record(
                    "CREST/xTB",
                    "running",
                    molecule_index=index,
                    molecule_total=len(seeds),
                    molecule_name=seed.molname,
                    generated_count=len(records),
                )
    else:
        if config.crest.enabled:
            warnings.append(
                "CREST/xTB requested but unavailable; using seed energies for smoke ranking."
            )
        records.extend(_crest_like_records_from_seeds(seeds, config))
        if progress is not None and seeds:
            progress.record(
                "CREST/xTB",
                "running",
                molecule_index=len(seeds),
                molecule_total=len(seeds),
                generated_count=len(records),
            )
    states.append(
        mark_done(
            step,
            seed_hash,
            config,
            details={"count": len(records), "enabled": config.crest.enabled},
        )
    )
    return records


def _crest_like_records_from_seeds(
    seeds: list[SeedConformerRecord],
    config: RunConfig,
) -> list[CrestConformerRecord]:
    records = []
    for index, seed in enumerate(seeds, start=1):
        metadata = {
            "ranking": {"source_workdir": None},
            "smoke_mode": True,
            "source_seed_id": seed.id,
        }
        records.append(
            CrestConformerRecord(
                id=f"{seed.id}_crest_smoke",
                parent_id=seed.id,
                input_molecule_id=seed.input_molecule_id,
                molname=seed.molname,
                canonical_smiles=seed.canonical_smiles,
                isomeric_smiles=seed.isomeric_smiles,
                molecular_formula=seed.molecular_formula,
                formal_charge=seed.formal_charge,
                explicit_proton_count=seed.explicit_proton_count,
                source_software="rdkit-seed-smoke-ranking",
                source_python_function="dsvr.workflow.engine._crest_like_records_from_seeds",
                warnings=["CREST/xTB skipped; RDKit seed energy used for smoke-mode ranking."],
                metadata=metadata,
                crest_index=index,
                energy_kcal_mol=seed.energy_kcal_mol if seed.energy_kcal_mol is not None else 0.0,
                relative_energy_kcal_mol=None,
            )
        )
    return records


def _run_optional_qm(
    ranked: list[RankedVariantRecord],
    config: RunConfig,
) -> list[RankedVariantRecord]:
    if not config.qm.enabled:
        return []
    ranked = sorted(ranked, key=lambda record: (record.rank, record.id))[: config.qm.max_candidates]
    backend = config.qm.backend if config.qm.backend != "none" else config.refinement.qm_backend
    if backend == "psi4" or config.refinement.psi4_enabled:
        return rescore_top_ranked_with_psi4(ranked, config)
    if backend == "pyscf" or config.refinement.pyscf_enabled:
        return rescore_top_ranked_with_pyscf(ranked, config)
    return []


def _select_thermo_inputs(
    records: list[CrestConformerRecord],
    config: RunConfig,
) -> list[CrestConformerRecord]:
    by_molecule: dict[str, list[CrestConformerRecord]] = {}
    for record in records:
        if record.energy_kcal_mol is None:
            continue
        by_molecule.setdefault(record.input_molecule_id, []).append(record)

    selected: list[CrestConformerRecord] = []
    for molecule_records in by_molecule.values():
        by_variant: dict[str | None, list[CrestConformerRecord]] = {}
        for record in sorted(molecule_records, key=lambda item: (item.energy_kcal_mol, item.id)):
            by_variant.setdefault(record.parent_id, []).append(record)
        ranked_variants = sorted(
            by_variant.items(),
            key=lambda item: (
                min(record.energy_kcal_mol or float("inf") for record in item[1]),
                str(item[0]),
            ),
        )[: config.thermo.max_variants_per_molecule]
        for _variant_id, variant_records in ranked_variants:
            selected.extend(
                sorted(
                    variant_records,
                    key=lambda item: (item.energy_kcal_mol or float("inf"), item.id),
                )[: config.thermo.max_conformers_per_variant]
            )
    return selected


def _load_protomers(outdir: Path) -> list[ProtomerRecord]:
    from dsvr.chemistry.tautomers import read_protomers_sdf

    records = []
    for path in sorted((outdir / "enumeration" / "protomers").glob("*_protomers.sdf")):
        records.extend(read_protomers_sdf(path))
    return records


def _load_tautomers(outdir: Path) -> list[Any]:
    from dsvr.chemistry.stereochemistry import read_tautomers_sdf

    records = []
    for path in sorted((outdir / "enumeration" / "tautomers").glob("*_tautomers.sdf")):
        records.extend(read_tautomers_sdf(path))
    return records


def _load_stereos(outdir: Path) -> list[Any]:
    from dsvr.chemistry.conformers_rdkit import read_stereo_sdf

    records = []
    for path in sorted((outdir / "enumeration" / "stereoisomers").glob("*_stereoisomers.sdf")):
        records.extend(read_stereo_sdf(path))
    return records


def _load_seeds(outdir: Path) -> list[SeedConformerRecord]:
    records: list[SeedConformerRecord] = []
    for path in sorted((outdir / "seeding").glob("**/*_seeds.sdf")):
        records.extend(read_seed_sdf(path))
    return records


def _load_crest_records(outdir: Path) -> list[CrestConformerRecord]:
    records: list[CrestConformerRecord] = []
    for path in sorted((outdir / "crest").glob("**/crest_provenance.jsonl")):
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.strip():
                continue
            try:
                records.append(CrestConformerRecord.model_validate(json.loads(line)))
            except (json.JSONDecodeError, ValueError):
                continue
    return records


def _load_thermo_records(outdir: Path) -> list[ThermoRecord]:
    records: list[ThermoRecord] = []
    for path in sorted((outdir / "xtb").glob("**/xtb_thermo.json")):
        try:
            records.append(ThermoRecord.model_validate_json(path.read_text(encoding="utf-8")))
        except ValueError:
            continue
    return records


def _load_existing_protomer_outputs(
    molecule: MoleculeInput,
    outdir: Path,
    config: RunConfig,
) -> list[ProtomerRecord]:
    if config.overwrite or not config.resume:
        return []
    path = outdir / "enumeration" / "protomers" / f"{molecule.input_id}_protomers.sdf"
    return _load_records_if_complete(path, read_protomers_sdf)


def _load_existing_tautomer_outputs(
    protomer: ProtomerRecord,
    outdir: Path,
    config: RunConfig,
) -> list[Any]:
    if config.overwrite or not config.resume:
        return []
    path = outdir / "enumeration" / "tautomers" / f"{protomer.id}_tautomers.sdf"
    return _load_records_if_complete(path, read_tautomers_sdf)


def _load_existing_stereo_outputs(
    tautomer: AnyLineageRecord,
    outdir: Path,
    config: RunConfig,
) -> list[Any]:
    if config.overwrite or not config.resume:
        return []
    path = outdir / "enumeration" / "stereoisomers" / f"{tautomer.id}_stereoisomers.sdf"
    return _load_records_if_complete(path, read_stereo_sdf)


def _load_existing_seed_outputs(
    stereo: AnyLineageRecord,
    outdir: Path,
    config: RunConfig,
) -> list[SeedConformerRecord]:
    if config.overwrite or not config.resume:
        return []
    path = outdir / "seeding" / "rdkit" / f"{stereo.id}_seeds.sdf"
    return _load_records_if_complete(path, read_seed_sdf)


def _load_records_if_complete(path: Path, loader: Callable[[Path], list[Any]]) -> list[Any]:
    companion_csv = path.with_suffix(".csv")
    if not path.exists() or not companion_csv.exists():
        return []
    try:
        records = loader(path)
    except (OSError, RuntimeError, ValueError):
        return []
    return records if records else []


def _write_input_table(path: Path, molecules: list[MoleculeInput]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "input_id",
                "molname",
                "source_format",
                "canonical_smiles",
                "isomeric_smiles",
            ],
        )
        writer.writeheader()
        for molecule in molecules:
            writer.writerow(
                {
                    "input_id": molecule.input_id,
                    "molname": molecule.molname,
                    "source_format": molecule.source_format,
                    "canonical_smiles": molecule.canonical_smiles,
                    "isomeric_smiles": molecule.isomeric_smiles,
                }
            )


def _ensure_invalid_inputs_csv(path: Path) -> None:
    if path.exists():
        return
    path.write_text("record_index,raw_record,error\n", encoding="utf-8")


def _write_protomer_outputs(record: ProtomerRecord) -> None:
    sdf_path, csv_path = record.output_paths
    writer = Chem.SDWriter(str(sdf_path))
    mol = Chem.Mol(record.rdkit_mol)
    mol.SetProp("_Name", record.id)
    mol.SetProp("DSVR_PROTOMER_ID", record.id)
    mol.SetProp("DSVR_INPUT_ID", record.input_molecule_id)
    mol.SetProp("DSVR_PARENT_ID", record.parent_id or "")
    mol.SetProp("DSVR_MOLNAME", record.molname)
    writer.write(mol)
    writer.close()
    with Path(csv_path).open("w", encoding="utf-8", newline="") as handle:
        writer_csv = csv.DictWriter(handle, fieldnames=["id", "parent_id", "warnings"])
        writer_csv.writeheader()
        writer_csv.writerow(
            {
                "id": record.id,
                "parent_id": record.parent_id,
                "warnings": " | ".join(record.warnings),
            }
        )


def _manifest(
    config: RunConfig,
    states: list[StepState],
    warnings: list[str],
    *,
    dry_run: bool,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "run_name": config.run_name,
        "input_path": str(config.input_path),
        "output_dir": str(config.output_dir),
        "dry_run": dry_run,
        "max_workers": config.max_workers,
        "tool_nproc": {"crest": config.crest.nproc},
        "steps": [
            {
                "name": state.step.name,
                "status": state.status,
                "skipped": state.skipped,
                "output_dir": str(state.output_dir),
                "details": state.details,
            }
            for state in states
        ],
        "warnings": warnings,
    } | (extra or {})


def _final_output_files(outdir: Path) -> list[Path]:
    names = [
        "manifest.json",
        "resolved_config.yaml",
        "invalid_inputs.csv",
        "inputs.csv",
        "protomers.csv",
        "tautomers.csv",
        "stereoisomers.csv",
        "seeds.csv",
        "crest_conformers.csv",
        "thermo.csv",
        "ranked_variants.csv",
        "ranked_variants.json",
        "ranked_variants.sdf",
        "report.md",
    ]
    return [outdir / name for name in names]


def _tool_available(executable: str) -> bool:
    return shutil.which(executable) is not None


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)
