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
from dsvr.chemistry.final3d import generate_final_3d_variants
from dsvr.chemistry.protonation import generate_protomer_candidates
from dsvr.chemistry.stereochemistry import enumerate_stereoisomers, read_tautomers_sdf
from dsvr.chemistry.stereo_auto3d_filter import filter_stereoisomers_with_auto3d
from dsvr.chemistry.tautomer_auto3d_filter import filter_tautomers_with_auto3d
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
from dsvr.workflow.recovery import (
    FailureKind,
    WorkflowRecoveryRecorder,
    classify_failure,
    safe_action_for_failure,
    should_skip_item_state,
)
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
    recovery = WorkflowRecoveryRecorder(outdir)
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
    for molecule in molecules:
        recovery.molecule(
            item_id=molecule.input_id,
            item_name=molecule.molname,
            stage="Input validation",
            status="completed",
        )
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
            if should_skip_item_state(
                recovery,
                molecule.input_id,
                resume=config.resume,
                stage="Protomer generation",
            ):
                progress.record(
                    "Protomer generation",
                    "running",
                    molecule_index=index,
                    molecule_total=len(molecules),
                    molecule_name=molecule.molname,
                    generated_count=len(protomers),
                    skipped_count=1,
                    message="Skipped previously checkpointed molecule.",
                )
                continue
            try:
                existing = _load_existing_protomer_outputs(molecule, outdir, config)
                if existing:
                    protomers.extend(existing)
                elif not config.protonation.enabled:
                    warning = "protonation.enabled=false; retained input state only"
                    warnings.append(warning)
                    protomers.extend(_fallback_protomer(molecule, config, warning=warning))
                else:
                    protomers.extend(generate_protomer_candidates(molecule, config))
                recovery.molecule(
                    item_id=molecule.input_id,
                    item_name=molecule.molname,
                    stage="Protomer generation",
                    status="completed",
                )
            except MolscrubUnavailableError as exc:
                _record_item_failure(
                    recovery,
                    progress,
                    config,
                    stage="Protomer generation",
                    item_id=molecule.input_id,
                    item_name=molecule.molname,
                    exc=exc,
                )
                if config.error_handling.keep_fallback_parent_state:
                    warning = f"Protomer generation failed; retained input state fallback: {exc}"
                    warnings.append(warning)
                    protomers.extend(_fallback_protomer(molecule, config, warning=warning))
                else:
                    recovery.molecule(
                        item_id=molecule.input_id,
                        item_name=molecule.molname,
                        stage="Protomer generation",
                        status="skipped",
                        message=str(exc),
                    )
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
    _write_stage_summary_sdf(outdir / "all_protomers.sdf", protomers)
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
    tautomer_input_hash = records_hash(protomers)
    if should_skip_step(steps["tautomers"], tautomer_input_hash, config):
        states.append(skipped_state(steps["tautomers"], tautomer_input_hash, config))
        tautomers = _load_tautomers(outdir)
        progress.record("Tautomer enumeration", "skipped", generated_count=len(tautomers))
    else:
        tautomers = []
        missing_protomers: list[ProtomerRecord] = []
        for protomer in protomers:
            existing = _load_existing_tautomer_outputs(protomer, outdir, config)
            if existing:
                tautomers.extend(existing)
            else:
                missing_protomers.append(protomer)
        if missing_protomers:
            if config.workflow_mode == "ligprep_like" and config.tautomer_filtering.enabled:
                try:
                    tautomers.extend(filter_tautomers_with_auto3d(missing_protomers, config))
                    for protomer in missing_protomers:
                        recovery.molecule(
                            item_id=protomer.id,
                            item_name=protomer.molname,
                            stage="Tautomer enumeration",
                            status="completed",
                        )
                except Exception as exc:
                    _record_item_failure(
                        recovery,
                        progress,
                        config,
                        stage="Tautomer enumeration",
                        item_id=None,
                        item_name=None,
                        exc=exc,
                    )
                    for protomer in missing_protomers:
                        if should_skip_item_state(
                            recovery,
                            protomer.id,
                            resume=config.resume,
                            stage="Tautomer enumeration",
                        ):
                            continue
                        try:
                            tautomers.extend(enumerate_tautomers(protomer, config))
                            recovery.molecule(
                                item_id=protomer.id,
                                item_name=protomer.molname,
                                stage="Tautomer enumeration",
                                status="completed",
                            )
                        except Exception as item_exc:
                            kind = classify_failure(item_exc, stage="Tautomer enumeration")
                            retry_config = _retry_reduced_enumeration_config(config, kind)
                            try:
                                tautomers.extend(enumerate_tautomers(protomer, retry_config))
                                recovery.molecule(
                                    item_id=protomer.id,
                                    item_name=protomer.molname,
                                    stage="Tautomer enumeration",
                                    status="completed",
                                    action=safe_action_for_failure(
                                        kind,
                                        stage="Tautomer enumeration",
                                    ),
                                    message=(
                                        f"Recovered after {kind.value} with "
                                        "reduced caps/timeouts."
                                    ),
                                )
                            except Exception as retry_exc:
                                _record_item_failure(
                                    recovery,
                                    progress,
                                    config,
                                    stage="Tautomer enumeration",
                                    item_id=protomer.id,
                                    item_name=protomer.molname,
                                    exc=retry_exc,
                                )
            else:
                for protomer in missing_protomers:
                    if should_skip_item_state(
                        recovery,
                        protomer.id,
                        resume=config.resume,
                        stage="Tautomer enumeration",
                    ):
                        continue
                    try:
                        tautomers.extend(enumerate_tautomers(protomer, config))
                        recovery.molecule(
                            item_id=protomer.id,
                            item_name=protomer.molname,
                            stage="Tautomer enumeration",
                            status="completed",
                        )
                    except Exception as exc:
                        kind = classify_failure(exc, stage="Tautomer enumeration")
                        retry_config = _retry_reduced_enumeration_config(config, kind)
                        try:
                            tautomers.extend(enumerate_tautomers(protomer, retry_config))
                            recovery.molecule(
                                item_id=protomer.id,
                                item_name=protomer.molname,
                                stage="Tautomer enumeration",
                                status="completed",
                                action=safe_action_for_failure(
                                    kind,
                                    stage="Tautomer enumeration",
                                ),
                                message=(
                                    f"Recovered after {kind.value} with reduced caps/timeouts."
                                ),
                            )
                        except Exception as retry_exc:
                            _record_item_failure(
                                recovery,
                                progress,
                                config,
                                stage="Tautomer enumeration",
                                item_id=protomer.id,
                                item_name=protomer.molname,
                                exc=retry_exc,
                            )
        states.append(
            mark_done(
                steps["tautomers"],
                tautomer_input_hash,
                config,
                details={"count": len(tautomers)},
            )
        )
    records.extend(tautomers)
    _write_stage_summary_sdf(outdir / "all_tautomers.sdf", tautomers)
    progress.record("Tautomer enumeration", "completed", generated_count=len(tautomers))

    progress.record("Stereoisomer enumeration", "started", generated_count=len(tautomers))
    stereos = _run_step_list(
        "stereochemistry",
        tautomers,
        config,
        states,
        lambda item: enumerate_stereoisomers(item, config),
        resume_loader=lambda: _load_stereos(outdir),
        retry_fn=lambda item, retry_config: enumerate_stereoisomers(item, retry_config),
        existing_loader=lambda item: _load_existing_stereo_outputs(item, outdir, config),
        progress=progress,
        progress_stage="Stereoisomer enumeration",
    )
    stereo_energy_result = filter_stereoisomers_with_auto3d(stereos, config)
    stereos = stereo_energy_result.all_records
    stereos_for_seeding = stereo_energy_result.selected_records
    records.extend(stereos)
    _write_stage_summary_sdf(outdir / "all_stereoisomers.sdf", stereos)
    progress.record(
        "Stereoisomer enumeration",
        "completed",
        generated_count=len(stereos),
        accepted_count=len(stereos_for_seeding),
        rejected_count=len(stereo_energy_result.rejected_records),
    )
    progress.record("Cheap variant scoring", "started", generated_count=len(stereos_for_seeding))
    stereos_for_seeding, decisions = select_stereo_records(stereos_for_seeding, config, "pre_3d")
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

    if _use_final_auto3d_default_path(config):
        return _run_final_auto3d_default_path(
            config=config,
            molecules=molecules,
            stereos_for_seeding=stereos_for_seeding,
            stereo_energy_result=stereo_energy_result,
            records=records,
            states=states,
            filtering_decisions=filtering_decisions,
            warnings=warnings,
            progress=progress,
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
                    active_command="Auto3D",
                )
                seeds.extend(
                    _generate_auto3d_seeds_with_retries(
                        stereos_for_seeding,
                        seed_config,
                        progress,
                    )
                )
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
    _write_stage_summary_sdf(outdir / "all_3d_conformers.sdf", seeds)
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

    progress.record("Report writing", "started", generated_count=len(ranked))
    write_all_provenance_outputs(records, outdir)
    write_audit_tables(
        outdir,
        records,
        filtering_decisions,
        xtb_prefilter_decisions,
        stereo_reduction,
    )
    write_final_ranked_outputs(outdir, ranked, config)
    _publish_top_level_run_outputs(outdir)
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
        "stereo_energy_filtering": {
            "enabled": config.stereoisomer_filtering.enabled,
            "enumerated_count": len(stereo_energy_result.all_records),
            "selected_count": len(stereo_energy_result.selected_records),
            "rejected_count": len(stereo_energy_result.rejected_records),
            "collapsed_count": stereo_energy_result.collapsed_count,
            "energy_evaluation_count": stereo_energy_result.energy_evaluation_count,
        },
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
    _publish_top_level_run_outputs(outdir)
    _record_progress_warnings(progress, warnings, "Report writing")
    progress.record("Report writing", "completed", generated_count=len(ranked))
    return WorkflowResult(
        outdir=outdir,
        molecule_count=len(molecules),
        ranked_records=legacy_records,
    )


def _record_progress_warnings(
    progress: ProgressRecorder,
    warnings: list[str],
    stage: str,
) -> None:
    seen: set[str] = set()
    for warning in warnings:
        if warning in seen:
            continue
        seen.add(warning)
        progress.warning(stage, warning)



def _record_item_failure(
    recovery: WorkflowRecoveryRecorder,
    progress: ProgressRecorder | None,
    config: RunConfig,
    *,
    stage: str,
    item_id: str | None,
    item_name: str | None,
    exc: BaseException,
) -> None:
    kind = classify_failure(exc, stage=stage)
    action = safe_action_for_failure(kind, stage=stage)
    recovery.failure(stage=stage, item_id=item_id, item_name=item_name, exc=exc, action=action)
    if progress is not None:
        progress.failure(stage, f"{kind.value}: {exc}")
    if config.error_handling.fail_fast:
        raise exc
    molecule_failure = kind in {
        FailureKind.INPUT_ERROR,
        FailureKind.PROTOMER_GENERATION_ERROR,
        FailureKind.DISK_LIMIT,
    }
    if molecule_failure and not config.error_handling.skip_failed_molecule:
        raise exc
    if not molecule_failure and not config.error_handling.skip_failed_variant:
        raise exc


def _retry_reduced_enumeration_config(config: RunConfig, kind: FailureKind) -> RunConfig:
    data = config.model_dump(mode="python")
    if (
        kind == FailureKind.TAUTOMER_TIMEOUT
        and config.error_handling.reduce_tautomer_cap_on_timeout
    ):
        data["enumeration"]["max_tautomers_per_protomer"] = max(
            1,
            config.enumeration.max_tautomers_per_protomer // 2,
        )
        data["enumeration"]["tautomer_timeout_seconds"] = max(
            1,
            config.enumeration.tautomer_timeout_seconds // 2,
        )
    if (
        kind == FailureKind.STEREO_TIMEOUT
        and config.error_handling.reduce_stereo_cap_on_timeout
    ):
        data["enumeration"]["max_stereoisomers_per_tautomer"] = max(
            1,
            config.enumeration.max_stereoisomers_per_tautomer // 2,
        )
        data["stereoisomer_filtering"]["timeout_seconds_per_tautomer"] = max(
            1,
            config.stereoisomer_filtering.timeout_seconds_per_tautomer // 2,
        )
    return RunConfig.model_validate(data)


def run_smoke_workflow(config: RunConfig) -> WorkflowResult:
    return run_workflow(config)


def _use_final_auto3d_default_path(config: RunConfig) -> bool:
    return (
        config.workflow_mode == "ligprep_like"
        and config.protocol == "default"
        and config.final_3d.tool == "auto3d"
        and config.final_3d.one_conformer_per_variant
    )


def _run_final_auto3d_default_path(
    *,
    config: RunConfig,
    molecules: list[MoleculeInput],
    stereos_for_seeding: list[AnyLineageRecord],
    stereo_energy_result: Any,
    records: list[AnyLineageRecord],
    states: list[StepState],
    filtering_decisions: list[FilteringDecision],
    warnings: list[str],
    progress: ProgressRecorder,
) -> WorkflowResult:
    outdir = config.output_dir
    steps = {step.name: step for step in planned_steps(config)}

    progress.record("Final Auto3D 3D generation", "started", generated_count=len(stereos_for_seeding))
    final_result = generate_final_3d_variants(stereos_for_seeding, config)
    warnings.extend(final_result.warnings)
    final_records = final_result.records
    states.append(
        mark_done(
            steps["seeding"],
            records_hash(stereos_for_seeding),
            config,
            details={
                "count": len(final_records),
                "method": "final_3d_auto3d",
                "one_conformer_per_variant": config.final_3d.one_conformer_per_variant,
                "used_fallback": final_result.used_fallback,
            },
        )
    )
    records.extend(final_records)
    _write_stage_summary_sdf(outdir / "all_3d_conformers.sdf", final_records)
    progress.record(
        "Final Auto3D 3D generation",
        "completed",
        generated_count=len(final_records),
        accepted_count=len(final_records),
    )

    stereo_reduction = StereoReductionResult(
        selected_seeds=final_records,
        decisions=[],
        representative_to_equivalent_seed_ids={},
        jobs_saved=0,
    )
    xtb_prefilter_decisions: list[XtbPrefilterDecision] = []
    _write_filtering_outputs(
        outdir,
        filtering_decisions,
        stereo_reduction,
        xtb_prefilter_decisions,
    )

    progress.record("Ranking", "started", generated_count=len(final_records))
    ranked = compute_delta_g_and_populations(final_records, config)
    write_ranked_outputs(ranked, outdir / "ranking")
    states.append(
        mark_done(
            steps["ranking"],
            records_hash(final_records),
            config,
            details={"count": len(ranked), "source": "final_3d_auto3d"},
        )
    )
    records.extend(ranked)
    progress.record("Ranking", "completed", generated_count=len(ranked))

    validation_result = _run_optional_crest_validation(
        final_records=final_records,
        filtering_decisions=filtering_decisions,
        config=config,
        warnings=warnings,
        progress=progress,
    )
    records.extend(validation_result["crest_records"])
    records.extend(validation_result["thermo_records"])
    states.append(
        mark_done(
            steps["crest"],
            records_hash(validation_result["selected_records"]),
            config,
            details={
                "count": validation_result["crest_record_count"],
                "selected_final_variant_count": validation_result["selected_count"],
                "enabled": config.optional_validation.crest_xtb_enabled,
                "optional_validation": True,
            },
        )
    )
    states.append(
        mark_done(
            steps["xtb_thermo"],
            records_hash(validation_result["crest_records"]),
            config,
            details={
                "count": validation_result["thermo_record_count"],
                "input_count": validation_result["thermo_input_count"],
                "enabled": config.optional_validation.xtb_thermo_enabled,
                "optional_validation": True,
            },
        )
    )
    states.append(
        mark_done(
            steps["censo"],
            records_hash(ranked),
            config,
            details={"count": 0, "enabled": False},
        )
    )
    states.append(
        mark_done(
            steps["qm"],
            records_hash(ranked),
            config,
            details={"count": 0, "backend": config.refinement.qm_backend},
        )
    )

    progress.record("Report writing", "started", generated_count=len(ranked))
    write_all_provenance_outputs(records, outdir)
    write_audit_tables(
        outdir,
        records,
        filtering_decisions,
        xtb_prefilter_decisions,
        stereo_reduction,
    )
    write_final_ranked_outputs(outdir, ranked, config)
    _publish_top_level_run_outputs(outdir)
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
        "stereo_energy_filtering": {
            "enabled": config.stereoisomer_filtering.enabled,
            "enumerated_count": len(stereo_energy_result.all_records),
            "selected_count": len(stereo_energy_result.selected_records),
            "rejected_count": len(stereo_energy_result.rejected_records),
            "collapsed_count": stereo_energy_result.collapsed_count,
            "energy_evaluation_count": stereo_energy_result.energy_evaluation_count,
        },
        "final_3d": {
            "enabled": True,
            "tool": config.final_3d.tool,
            "requested_count": len(stereos_for_seeding),
            "final_conformer_count": len(final_records),
            "one_conformer_per_variant": config.final_3d.one_conformer_per_variant,
            "used_fallback": final_result.used_fallback,
        },
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
    manifest["optional_validation"] = {
        "crest_xtb_enabled": config.optional_validation.crest_xtb_enabled,
        "xtb_thermo_enabled": config.optional_validation.xtb_thermo_enabled,
        "selection": config.optional_validation.selection,
        "max_variants_per_molecule": config.optional_validation.max_variants_per_molecule,
        "selected_count": validation_result["selected_count"],
        "crest_record_count": validation_result["crest_record_count"],
        "thermo_record_count": validation_result["thermo_record_count"],
        "outputs": validation_result["outputs"],
        "ranking_overwritten": False,
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
    _publish_top_level_run_outputs(outdir)
    _record_progress_warnings(progress, warnings, "Report writing")
    progress.record("Report writing", "completed", generated_count=len(ranked))
    return WorkflowResult(
        outdir=outdir,
        molecule_count=len(molecules),
        ranked_records=legacy_records,
    )


def _progress_stage_names(config: RunConfig) -> list[str]:
    if config.protocol == "auto3d_entropy":
        return [
            "Input validation",
            "Standardization",
            "Protomer generation",
            "Auto3D representative generation",
            "Representative plausibility scoring",
            "Ranking",
            "Report writing",
        ]
    if _use_final_auto3d_default_path(config):
        stages = [
            "Input validation",
            "Standardization",
            "Protomer generation",
            "Tautomer enumeration",
            "Stereoisomer enumeration",
            "Cheap variant scoring",
            "Final Auto3D 3D generation",
            "Ranking",
        ]
        if config.optional_validation.crest_xtb_enabled or config.optional_validation.xtb_thermo_enabled:
            stages.append("Optional validation")
        stages.append("Report writing")
        return stages
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
        "Report writing",
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

    _write_stage_summary_sdf(outdir / "all_tautomers.sdf", [])
    _write_stage_summary_sdf(outdir / "all_stereoisomers.sdf", [])

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
    _write_stage_summary_sdf(outdir / "all_3d_conformers.sdf", seeds)
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

    progress.record("Report writing", "started", generated_count=len(ranked))
    write_all_provenance_outputs(records, outdir)
    write_audit_tables(outdir, records, [], [], empty_stereo_reduction)
    write_final_ranked_outputs(outdir, ranked, config)
    _write_auto3d_protocol_structure_summary(
        outdir,
        molecules=molecules,
        protomers=protomers,
        seeds=seeds,
        score_records=score_records,
        ranked=ranked,
    )
    _publish_top_level_run_outputs(outdir)
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
    _publish_top_level_run_outputs(outdir)
    _record_progress_warnings(progress, warnings, "Report writing")
    progress.record("Report writing", "completed", generated_count=len(ranked))
    return WorkflowResult(
        outdir=outdir,
        molecule_count=len(molecules),
        ranked_records=legacy_records,
    )


def _write_auto3d_protocol_structure_summary(
    outdir: Path,
    *,
    molecules: list[MoleculeInput],
    protomers: list[ProtomerRecord],
    seeds: list[SeedConformerRecord],
    score_records: list[ThermoRecord],
    ranked: list[RankedVariantRecord],
) -> None:
    seed_by_protomer: dict[str, list[SeedConformerRecord]] = {}
    fallback_seeds = []
    for seed in seeds:
        if seed.parent_id:
            seed_by_protomer.setdefault(seed.parent_id, []).append(seed)
        if "auto3d_fallback" in seed.metadata:
            fallback_seeds.append(seed)

    missing_protomers = [protomer for protomer in protomers if protomer.id not in seed_by_protomer]
    molecule_ids_with_seeds = {seed.input_molecule_id for seed in seeds}
    molecules_without_structures = [
        molecule for molecule in molecules if molecule.input_id not in molecule_ids_with_seeds
    ]
    failure_rows = _auto3d_structure_failure_rows(fallback_seeds, missing_protomers)

    summary_rows = [
        {
            "stage": "Input validation",
            "generated_structures": len(molecules),
            "input_molecules": len(molecules),
            "molecules_with_structures": len(molecules),
            "failed_molecules": 0,
            "fallback_structures": 0,
            "output_sdf": "",
            "output_csv": "input/inputs.csv",
            "notes": "validated input records",
        },
        {
            "stage": "Protomer generation",
            "generated_structures": len(protomers),
            "input_molecules": len(molecules),
            "molecules_with_structures": len({item.input_molecule_id for item in protomers}),
            "failed_molecules": len(molecules)
            - len({item.input_molecule_id for item in protomers}),
            "fallback_structures": 0,
            "output_sdf": "all_protomers.sdf",
            "output_csv": "enumeration/protomers/*_protomers.csv",
            "notes": (
                "molscrub protomer candidates, or original molecule fallback "
                "if molscrub unavailable"
            ),
        },
        {
            "stage": "Auto3D representative generation",
            "generated_structures": len(seeds),
            "input_molecules": len(molecules),
            "molecules_with_structures": len(molecule_ids_with_seeds),
            "failed_molecules": len(molecules_without_structures),
            "fallback_structures": len(fallback_seeds),
            "output_sdf": "all_3d_conformers.sdf",
            "output_csv": "auto3d_protocol_seeds.csv",
            "notes": (
                f"{len(fallback_seeds)} protomer(s) used RDKit fallback after Auto3D failed; "
                f"{len(missing_protomers)} protomer(s) have no generated structure"
            ),
        },
        {
            "stage": "Representative plausibility scoring",
            "generated_structures": len(score_records),
            "input_molecules": len(molecules),
            "molecules_with_structures": len({item.input_molecule_id for item in score_records}),
            "failed_molecules": len(molecules)
            - len({item.input_molecule_id for item in score_records}),
            "fallback_structures": 0,
            "output_sdf": "",
            "output_csv": "auto3d_representative_scores.csv",
            "notes": "SVPScore-style representative scoring records",
        },
        {
            "stage": "Ranking",
            "generated_structures": len(ranked),
            "input_molecules": len(molecules),
            "molecules_with_structures": len({item.input_molecule_id for item in ranked}),
            "failed_molecules": len(molecules) - len({item.input_molecule_id for item in ranked}),
            "fallback_structures": 0,
            "output_sdf": "ranked_variants.sdf",
            "output_csv": "ranked_variants.csv",
            "notes": "final ranked structural variants",
        },
    ]
    _write_csv(outdir / "structure_generation_summary.csv", summary_rows)
    _write_csv(outdir / "structure_failures.csv", failure_rows, empty_columns=[
        "input_molecule_id",
        "molname",
        "protomer_id",
        "stage",
        "failure_kind",
        "fallback_structure_generated",
        "reason",
    ])


def _auto3d_structure_failure_rows(
    fallback_seeds: list[SeedConformerRecord],
    missing_protomers: list[ProtomerRecord],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for seed in fallback_seeds:
        fallback = seed.metadata.get("auto3d_fallback", {})
        rows.append(
            {
                "input_molecule_id": seed.input_molecule_id,
                "molname": seed.molname,
                "protomer_id": seed.parent_id,
                "stage": "Auto3D representative generation",
                "failure_kind": "auto3d_optimized_structure_failed",
                "fallback_structure_generated": True,
                "reason": fallback.get("reason") or "Auto3D produced no optimized representative",
            }
        )
    for protomer in missing_protomers:
        rows.append(
            {
                "input_molecule_id": protomer.input_molecule_id,
                "molname": protomer.molname,
                "protomer_id": protomer.id,
                "stage": "Auto3D representative generation",
                "failure_kind": "no_structure_generated",
                "fallback_structure_generated": False,
                "reason": "Neither Auto3D nor RDKit fallback produced a structure",
            }
        )
    return rows


def _publish_top_level_run_outputs(outdir: Path) -> None:
    link_sources = [
        outdir / "seeding" / "auto3d_protocol" / "auto3d_protocol_seeds.sdf",
        outdir / "seeding" / "auto3d_protocol" / "auto3d_protocol_seeds.csv",
        outdir / "seeding" / "auto3d_protocol" / "auto3d_adaptive_plan.csv",
        outdir / "auto3d_representatives" / "auto3d_representative_scores.csv",
        outdir / "ranking" / "ranking_summary.md",
    ]
    for source in link_sources:
        _ensure_top_level_link(outdir, source)

    inventory_targets = [
        outdir / "stage_summary.csv",
        outdir / "structure_generation_summary.csv",
        outdir / "structure_failures.csv",
        outdir / "ranked_variants.sdf",
        outdir / "ranked_variants.csv",
        outdir / "ranked_variants.json",
        outdir / "final_variants.sdf",
        outdir / "final_variants.csv",
        outdir / "final_variants.json",
        outdir / "final_variant_energies.csv",
        outdir / "all_protomers.sdf",
        outdir / "all_tautomers.sdf",
        outdir / "all_stereoisomers.sdf",
        outdir / "all_3d_conformers.sdf",
        outdir / "stereoisomers_all.csv",
        outdir / "stereoisomers_selected.csv",
        outdir / "stereoisomers_rejected.csv",
        outdir / "stereo_energy_ranked.csv",
        outdir / "stereo_enantiomer_groups.csv",
        outdir / "auto3d_protocol_seeds.sdf",
        outdir / "auto3d_protocol_seeds.csv",
        outdir / "auto3d_adaptive_plan.csv",
        outdir / "auto3d_representative_scores.csv",
        outdir / "ranking_summary.md",
        outdir / "summary.md",
        outdir / "report.md",
        outdir / "manifest.json",
    ]
    rows = []
    for target in inventory_targets:
        rows.append(
            {
                "artifact": target.name,
                "path": str(target.relative_to(outdir)),
                "kind": _artifact_kind(target),
                "exists": target.exists(),
                "size_bytes": target.stat().st_size if target.exists() else 0,
                "record_count": _artifact_record_count(target),
                "target": str(target.resolve().relative_to(outdir.resolve()))
                if target.exists() and target.is_symlink()
                else "",
            }
        )
    _write_csv(outdir / "run_outputs.csv", rows)


def _ensure_top_level_link(outdir: Path, source: Path) -> None:
    if not source.exists() or source.parent == outdir:
        return
    destination = outdir / source.name
    if destination.exists() or destination.is_symlink():
        if destination.is_dir() and not destination.is_symlink():
            return
        destination.unlink()
    try:
        destination.symlink_to(source.relative_to(outdir))
    except OSError:
        shutil.copy2(source, destination)


def _artifact_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".sdf":
        return "structure_sdf"
    if suffix == ".csv":
        return "table_csv"
    if suffix == ".json":
        return "metadata_json"
    if suffix == ".md":
        return "report_markdown"
    return "artifact"


def _artifact_record_count(path: Path) -> int | str:
    if not path.exists():
        return 0
    if path.suffix.lower() == ".csv":
        try:
            with path.open(encoding="utf-8") as handle:
                return max(0, sum(1 for _ in handle) - 1)
        except OSError:
            return "unknown"
    if path.suffix.lower() == ".sdf":
        return _count_sdf_records(path)
    return ""


def _count_sdf_records(path: Path) -> int | str:
    try:
        supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
        return sum(1 for molecule in supplier if molecule is not None)
    except (OSError, RuntimeError):
        return "unknown"


def _write_csv(
    path: Path,
    rows: list[dict[str, object]],
    *,
    empty_columns: list[str] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(rows[0].keys()) if rows else empty_columns or []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _run_step_list(
    step_name: str,
    inputs: list[Any],
    config: RunConfig,
    states: list[StepState],
    fn,
    resume_loader=None,
    existing_loader: Callable[[Any], list[Any]] | None = None,
    retry_fn: Callable[[Any, RunConfig], list[Any]] | None = None,
    progress: ProgressRecorder | None = None,
    progress_stage: str | None = None,
) -> list[Any]:
    step = {item.name: item for item in planned_steps(config)}[step_name]
    stage = progress_stage or step.description
    recovery = WorkflowRecoveryRecorder(config.output_dir)
    input_hash = records_hash(inputs)
    if should_skip_step(step, input_hash, config):
        states.append(skipped_state(step, input_hash, config))
        recovery.stage(stage, "skipped")
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
        item_id = getattr(item, "id", None) or getattr(item, "input_id", None) or str(index)
        item_name = getattr(item, "molname", None) or str(item_id)
        if should_skip_item_state(recovery, str(item_id), resume=config.resume, stage=stage):
            if progress is not None and progress_stage is not None:
                progress.record(
                    progress_stage,
                    "running",
                    molecule_index=index,
                    molecule_total=len(inputs),
                    molecule_name=item_name,
                    generated_count=len(outputs),
                    skipped_count=1,
                    message="Skipped previously checkpointed item.",
                )
            continue
        try:
            existing = existing_loader(item) if existing_loader is not None else []
            if existing:
                outputs.extend(existing)
            else:
                outputs.extend(fn(item))
            recovery.molecule(
                item_id=str(item_id),
                item_name=item_name,
                stage=stage,
                status="completed",
            )
        except Exception as exc:
            kind = classify_failure(exc, stage=stage)
            if (
                kind in {FailureKind.TAUTOMER_TIMEOUT, FailureKind.STEREO_TIMEOUT}
                and retry_fn is not None
            ):
                try:
                    outputs.extend(
                        retry_fn(item, _retry_reduced_enumeration_config(config, kind))
                    )
                    recovery.molecule(
                        item_id=str(item_id),
                        item_name=item_name,
                        stage=stage,
                        status="completed",
                        action=safe_action_for_failure(kind, stage=stage),
                        message=f"Recovered after {kind.value} with reduced caps/timeouts.",
                    )
                except Exception as retry_exc:
                    _record_item_failure(
                        recovery,
                        progress,
                        config,
                        stage=stage,
                        item_id=str(item_id),
                        item_name=item_name,
                        exc=retry_exc,
                    )
            else:
                _record_item_failure(
                    recovery,
                    progress,
                    config,
                    stage=stage,
                    item_id=str(item_id),
                    item_name=item_name,
                    exc=exc,
                )
        if progress is not None and progress_stage is not None:
            progress.record(
                progress_stage,
                "running",
                molecule_index=index,
                molecule_total=len(inputs),
                molecule_name=item_name,
                generated_count=len(outputs),
            )
    states.append(mark_done(step, input_hash, config, details={"count": len(outputs)}))
    return outputs


def _generate_auto3d_seeds_with_retries(
    stereos: list[AnyLineageRecord],
    config: RunConfig,
    progress: ProgressRecorder | None = None,
) -> list[SeedConformerRecord]:
    recovery = WorkflowRecoveryRecorder(config.output_dir)
    try:
        seeds = generate_auto3d_seeds(stereos, config)
        for stereo in stereos:
            recovery.molecule(
                item_id=stereo.id,
                item_name=stereo.molname,
                stage="3D seeding",
                status="completed",
            )
        return seeds
    except (Auto3DExecutionError, Auto3DUnavailableError) as exc:
        if config.seeding.auto3d_use_gpu and config.error_handling.retry_auto3d_cpu_on_gpu_failure:
            retry_data = config.model_dump(mode="python")
            retry_data["seeding"]["auto3d_use_gpu"] = False
            retry_config = RunConfig.model_validate(retry_data)
            try:
                seeds = generate_auto3d_seeds(stereos, retry_config)
                for stereo in stereos:
                    recovery.molecule(
                        item_id=stereo.id,
                        item_name=stereo.molname,
                        stage="3D seeding",
                        status="completed",
                        action="retry_auto3d_cpu_after_gpu_failure",
                        message=str(exc),
                    )
                return seeds
            except (Auto3DExecutionError, Auto3DUnavailableError) as cpu_exc:
                exc = cpu_exc
        if (
            config.error_handling.retry_auto3d_smaller_batch_on_batch_failure
            and len(stereos) > 1
        ):
            recovered: list[SeedConformerRecord] = []
            for index, stereo in enumerate(stereos, start=1):
                try:
                    recovered.extend(generate_auto3d_seeds([stereo], config))
                    recovery.molecule(
                        item_id=stereo.id,
                        item_name=stereo.molname,
                        stage="3D seeding",
                        status="completed",
                        action="retry_auto3d_smaller_batch_after_batch_failure",
                    )
                except (Auto3DExecutionError, Auto3DUnavailableError) as item_exc:
                    _record_item_failure(
                        recovery,
                        progress,
                        config,
                        stage="3D seeding",
                        item_id=stereo.id,
                        item_name=stereo.molname,
                        exc=item_exc,
                    )
                if progress is not None:
                    progress.record(
                        "3D seeding",
                        "running",
                        molecule_index=index,
                        molecule_total=len(stereos),
                        molecule_name=stereo.molname,
                        generated_count=len(recovered),
                        active_command="Auto3D",
                    )
            if recovered:
                return recovered
        raise exc


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


def _fallback_protomer(
    molecule: MoleculeInput,
    config: RunConfig,
    *,
    warning: str = "molscrub unavailable; input protomer fallback retained",
) -> list[ProtomerRecord]:
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


def _run_optional_crest_validation(
    *,
    final_records: list[SeedConformerRecord],
    filtering_decisions: list[FilteringDecision],
    config: RunConfig,
    warnings: list[str],
    progress: ProgressRecorder,
) -> dict[str, Any]:
    validation_dir = config.output_dir / "optional_validation"
    outputs = [
        config.output_dir / "crest_validation.csv",
        config.output_dir / "crest_validation.sdf",
        config.output_dir / "crest_validation_report.md",
    ]
    if not config.optional_validation.crest_xtb_enabled:
        return {
            "selected_records": [],
            "selected_count": 0,
            "crest_records": [],
            "crest_record_count": 0,
            "thermo_records": [],
            "thermo_record_count": 0,
            "thermo_input_count": 0,
            "outputs": [],
        }

    selected = _select_optional_validation_records(final_records, filtering_decisions, config)
    progress.record("Optional validation", "started", generated_count=len(selected))
    validation_dir.mkdir(parents=True, exist_ok=True)
    _write_optional_validation_input_sdf(validation_dir / "selected_final_variants.sdf", selected)
    _write_optional_validation_input_csv(validation_dir / "selected_final_variants.csv", selected)

    validation_config = _config_for_optional_validation(config)
    crest_records: list[CrestConformerRecord] = []
    crest_tools_available = _tool_available(config.crest.executable) and _tool_available(
        config.crest.xtb_executable
    )
    if selected and crest_tools_available and config.crest.enabled:
        for index, record in enumerate(selected, start=1):
            crest_records.extend(run_crest_for_seed(record, validation_config))
            progress.record(
                "Optional validation",
                "running",
                molecule_index=index,
                molecule_total=len(selected),
                molecule_name=record.molname,
                generated_count=len(crest_records),
                active_command="crest/xtb",
            )
    elif selected:
        message = "Optional CREST/xTB validation requested but CREST/xTB is unavailable."
        warnings.append(message)
        crest_records = _crest_like_records_from_seeds(selected, validation_config)

    thermo_records: list[ThermoRecord] = []
    thermo_inputs: list[CrestConformerRecord] = []
    if config.optional_validation.xtb_thermo_enabled and crest_records:
        thermo_inputs = _select_thermo_inputs(crest_records, validation_config)
        if _tool_available(config.crest.xtb_executable):
            for conformer in thermo_inputs:
                thermo_records.append(run_xtb_thermo(conformer, validation_config))
        else:
            warnings.append("Optional xTB thermo validation requested but xTB is unavailable.")

    _write_crest_validation_csv(outputs[0], selected, crest_records, thermo_records)
    _write_crest_validation_sdf(outputs[1], selected, crest_records)
    _write_crest_validation_report(outputs[2], config, selected, crest_records, thermo_records)
    progress.record(
        "Optional validation",
        "completed",
        generated_count=len(crest_records),
        accepted_count=len(selected),
    )
    return {
        "selected_records": selected,
        "selected_count": len(selected),
        "crest_records": crest_records,
        "crest_record_count": len(crest_records),
        "thermo_records": thermo_records,
        "thermo_record_count": len(thermo_records),
        "thermo_input_count": len(thermo_inputs),
        "outputs": [str(path) for path in outputs],
    }


def _select_optional_validation_records(
    final_records: list[SeedConformerRecord],
    filtering_decisions: list[FilteringDecision],
    config: RunConfig,
) -> list[SeedConformerRecord]:
    rescue_parent_ids = {
        decision.record_id
        for decision in filtering_decisions
        if decision.selected and decision.rescue_reason
    }
    grouped: dict[str, list[SeedConformerRecord]] = {}
    for record in final_records:
        grouped.setdefault(record.input_molecule_id, []).append(record)

    selected: list[SeedConformerRecord] = []
    for input_id in sorted(grouped):
        records = sorted(grouped[input_id], key=lambda item: _seed_energy_sort_key(item))
        chosen: list[SeedConformerRecord] = records[:3]
        for record in records:
            if record.parent_id in rescue_parent_ids and record not in chosen:
                chosen.append(record)
            if len(chosen) >= config.optional_validation.max_variants_per_molecule:
                break
        selected.extend(chosen[: config.optional_validation.max_variants_per_molecule])
    return selected


def _seed_energy_sort_key(record: SeedConformerRecord) -> tuple[float, str]:
    return (float("inf") if record.energy_kcal_mol is None else record.energy_kcal_mol, record.id)


def _config_for_optional_validation(config: RunConfig) -> RunConfig:
    data = config.model_dump(mode="python")
    keep_raw = config.optional_validation.keep_raw_xyz or config.optional_validation.cleanup_policy == "debug_all"
    compact = config.optional_validation.cleanup_policy == "compact"
    data["crest"]["keep_raw_xyz"] = keep_raw
    data["crest"]["compress_raw_outputs"] = compact or data["crest"].get("compress_raw_outputs", True)
    data["crest"]["delete_intermediate_xyz"] = compact or data["crest"].get("delete_intermediate_xyz", True)
    data["disk"]["keep_raw_xyz"] = keep_raw
    data["disk"]["compress_raw_outputs"] = compact or data["disk"].get("compress_raw_outputs", True)
    data["disk"]["delete_intermediate_xyz"] = compact or data["disk"].get("delete_intermediate_xyz", True)
    if config.optional_validation.xtb_thermo_enabled:
        data["thermo"]["enabled"] = True
        data["thermo"]["xtb_hessian"] = True
        data["thermo"]["xtb_thermo"] = True
    return RunConfig.model_validate(data)


def _write_optional_validation_input_sdf(path: Path, records: list[SeedConformerRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        if record.rdkit_mol is None:
            continue
        mol = Chem.Mol(record.rdkit_mol)
        mol.SetProp("_Name", record.id)
        mol.SetProp("DSVR_FINAL_VARIANT_ID", record.id)
        mol.SetProp("DSVR_STEREO_ID", record.parent_id or "")
        mol.SetProp("DSVR_INPUT_ID", record.input_molecule_id)
        mol.SetProp("DSVR_OPTIONAL_VALIDATION_SELECTED", "True")
        writer.write(mol)
    writer.close()


def _write_optional_validation_input_csv(path: Path, records: list[SeedConformerRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "final_variant_id",
                "input_id",
                "molname",
                "stereo_id",
                "auto3d_energy_kcal_mol",
            ],
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "final_variant_id": record.id,
                    "input_id": record.input_molecule_id,
                    "molname": record.molname,
                    "stereo_id": record.parent_id,
                    "auto3d_energy_kcal_mol": record.energy_kcal_mol,
                }
            )


def _write_crest_validation_csv(
    path: Path,
    selected: list[SeedConformerRecord],
    crest_records: list[CrestConformerRecord],
    thermo_records: list[ThermoRecord],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    selected_by_id = {record.id: record for record in selected}
    thermo_by_parent = {record.parent_id: record for record in thermo_records}
    columns = [
        "optional_validation",
        "input_id",
        "molname",
        "final_variant_id",
        "stereo_id",
        "auto3d_energy_kcal_mol",
        "crest_conformer_id",
        "crest_index",
        "crest_energy_kcal_mol",
        "crest_relative_energy_kcal_mol",
        "xtb_free_energy_kcal_mol",
        "xtb_entropy_cal_mol_k",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        if crest_records:
            for record in crest_records:
                seed = selected_by_id.get(record.parent_id or "")
                thermo = thermo_by_parent.get(record.id)
                writer.writerow(
                    {
                        "optional_validation": True,
                        "input_id": record.input_molecule_id,
                        "molname": record.molname,
                        "final_variant_id": record.parent_id,
                        "stereo_id": seed.parent_id if seed else None,
                        "auto3d_energy_kcal_mol": seed.energy_kcal_mol if seed else None,
                        "crest_conformer_id": record.id,
                        "crest_index": record.crest_index,
                        "crest_energy_kcal_mol": record.energy_kcal_mol,
                        "crest_relative_energy_kcal_mol": record.relative_energy_kcal_mol,
                        "xtb_free_energy_kcal_mol": thermo.free_energy_kcal_mol if thermo else None,
                        "xtb_entropy_cal_mol_k": thermo.entropy_cal_mol_k if thermo else None,
                        "warnings": " | ".join(record.warnings + (thermo.warnings if thermo else [])),
                    }
                )
        else:
            for record in selected:
                writer.writerow(
                    {
                        "optional_validation": True,
                        "input_id": record.input_molecule_id,
                        "molname": record.molname,
                        "final_variant_id": record.id,
                        "stereo_id": record.parent_id,
                        "auto3d_energy_kcal_mol": record.energy_kcal_mol,
                        "warnings": "CREST/xTB validation produced no conformer records.",
                    }
                )


def _write_crest_validation_sdf(
    path: Path,
    selected: list[SeedConformerRecord],
    crest_records: list[CrestConformerRecord],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    by_seed: dict[str, list[CrestConformerRecord]] = {}
    for record in crest_records:
        by_seed.setdefault(record.parent_id or "", []).append(record)
    writer = Chem.SDWriter(str(path))
    for seed in selected:
        if seed.rdkit_mol is None:
            continue
        validation_records = sorted(by_seed.get(seed.id, []), key=lambda item: (item.crest_index, item.id))
        best = min(validation_records, key=lambda item: (item.energy_kcal_mol or float("inf"), item.id), default=None)
        mol = Chem.Mol(seed.rdkit_mol)
        mol.SetProp("_Name", seed.id)
        mol.SetProp("DSVR_OPTIONAL_VALIDATION", "CREST/xTB")
        mol.SetProp("DSVR_OPTIONAL_VALIDATION_DOES_NOT_SET_RANKING", "True")
        mol.SetProp("DSVR_FINAL_VARIANT_ID", seed.id)
        mol.SetProp("DSVR_STEREO_ID", seed.parent_id or "")
        mol.SetProp("DSVR_INPUT_ID", seed.input_molecule_id)
        mol.SetProp("DSVR_FINAL_AUTO3D_ENERGY_KCAL_MOL", "" if seed.energy_kcal_mol is None else str(seed.energy_kcal_mol))
        mol.SetProp("DSVR_CREST_CONFORMER_COUNT", str(len(validation_records)))
        if best is not None:
            mol.SetProp("DSVR_BEST_CREST_CONFORMER_ID", best.id)
            mol.SetProp("DSVR_BEST_CREST_ENERGY_KCAL_MOL", "" if best.energy_kcal_mol is None else str(best.energy_kcal_mol))
            mol.SetProp("DSVR_BEST_CREST_RELATIVE_ENERGY_KCAL_MOL", "" if best.relative_energy_kcal_mol is None else str(best.relative_energy_kcal_mol))
        writer.write(mol)
    writer.close()


def _write_crest_validation_report(
    path: Path,
    config: RunConfig,
    selected: list[SeedConformerRecord],
    crest_records: list[CrestConformerRecord],
    thermo_records: list[ThermoRecord],
) -> None:
    by_input: dict[str, int] = {}
    for record in selected:
        by_input[record.input_molecule_id] = by_input.get(record.input_molecule_id, 0) + 1
    lines = [
        "# Optional CREST/xTB Validation Report",
        "",
        "This validation is optional and does not overwrite the default ligand-prep ranking.",
        "",
        f"- Enabled: {config.optional_validation.crest_xtb_enabled}",
        f"- Selection: {config.optional_validation.selection}",
        f"- Max variants per molecule: {config.optional_validation.max_variants_per_molecule}",
        f"- Selected final variants: {len(selected)}",
        f"- CREST conformer records: {len(crest_records)}",
        f"- xTB thermo records: {len(thermo_records)}",
        f"- Cleanup policy: {config.optional_validation.cleanup_policy}",
        f"- Keep raw XYZ: {config.optional_validation.keep_raw_xyz}",
        "",
        "## Selected Variants Per Input",
        "",
        "| Input ID | Selected variants |",
        "| --- | ---: |",
        *[f"| {input_id} | {count} |" for input_id, count in sorted(by_input.items())],
        "",
        "## Outputs",
        "",
        "- crest_validation.csv",
        "- crest_validation.sdf",
        "- crest_validation_report.md",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
                    active_command="crest/xtb",
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


def _write_stage_summary_sdf(path: Path, records: list[AnyLineageRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    for record in records:
        molecule = getattr(record, "rdkit_mol", None)
        if molecule is None:
            continue
        mol = Chem.Mol(molecule)
        mol.SetProp("_Name", record.id)
        props = {
            "DSVR_STAGE": record.stage_name,
            "DSVR_RECORD_ID": record.id,
            "DSVR_PARENT_ID": record.parent_id or "",
            "DSVR_INPUT_ID": record.input_molecule_id,
            "DSVR_MOLNAME": record.molname,
            "DSVR_CANONICAL_SMILES": record.canonical_smiles or "",
            "DSVR_ISOMERIC_SMILES": record.isomeric_smiles or "",
            "DSVR_FORMULA": record.molecular_formula or "",
            "DSVR_FORMAL_CHARGE": "" if record.formal_charge is None else str(record.formal_charge),
            "DSVR_EXPLICIT_PROTON_COUNT": ""
            if record.explicit_proton_count is None
            else str(record.explicit_proton_count),
        }
        if isinstance(record, ProtomerRecord):
            props["DSVR_PROTOMER_ID"] = record.id
        elif record.stage_name == "tautomer":
            props["DSVR_TAUTOMER_ID"] = record.id
            props["DSVR_PARENT_PROTOMER_ID"] = record.parent_id or ""
        elif record.stage_name == "stereo":
            props["DSVR_STEREO_ID"] = record.id
            props["DSVR_PARENT_TAUTOMER_ID"] = record.parent_id or ""
        elif isinstance(record, SeedConformerRecord):
            props["DSVR_SEED_ID"] = record.id
            props["DSVR_PARENT_STEREO_ID"] = record.parent_id or ""
            props["DSVR_FORCEFIELD"] = record.forcefield or ""
            props["DSVR_FORCEFIELD_STATUS"] = record.forcefield_status
            props["DSVR_EMBEDDING_STATUS"] = record.embedding_status
            props["DSVR_ENERGY_KCAL_MOL"] = (
                "" if record.energy_kcal_mol is None else str(record.energy_kcal_mol)
            )
            props["DSVR_AUTO3D_FALLBACK"] = str("auto3d_fallback" in record.metadata)
        for key, value in props.items():
            mol.SetProp(key, value)
        writer.write(mol)
    writer.close()


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
        "stereoisomers_all.csv",
        "stereoisomers_selected.csv",
        "stereoisomers_rejected.csv",
        "stereo_energy_ranked.csv",
        "stereo_enantiomer_groups.csv",
        "all_protomers.sdf",
        "all_tautomers.sdf",
        "all_stereoisomers.sdf",
        "all_3d_conformers.sdf",
        "seeds.csv",
        "crest_conformers.csv",
        "thermo.csv",
        "ranked_variants.csv",
        "ranked_variants.json",
        "ranked_variants.sdf",
        "final_variants.sdf",
        "final_variants.csv",
        "final_variants.json",
        "final_variant_energies.csv",
        "variant_decisions.csv",
        "protomers_all.csv",
        "protomers_selected.csv",
        "protomers_rejected.csv",
        "tautomers_all_pre_auto3d.csv",
        "tautomers_auto3d_ranked.csv",
        "tautomers_selected.csv",
        "tautomers_rejected.csv",
        "crest_validation.csv",
        "crest_validation.sdf",
        "crest_validation_report.md",
        "report.md",
    ]
    return [outdir / name for name in names]


def _tool_available(executable: str) -> bool:
    return shutil.which(executable) is not None


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)
