from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

from rdkit.Chem import Lipinski

from dsvr.filtering.selection import FilteringDecision
from dsvr.filtering.stereo_reduce import StereoReductionResult
from dsvr.filtering.xtb_prefilter import XtbPrefilterDecision
from dsvr.models import AnyLineageRecord


def write_audit_tables(
    outdir: Path,
    records: list[AnyLineageRecord],
    filtering_decisions: list[FilteringDecision],
    xtb_prefilter_decisions: list[XtbPrefilterDecision],
    stereo_reduction: StereoReductionResult,
) -> None:
    write_enumeration_counts(outdir / "enumeration_counts.csv", records)
    write_variant_selection(
        outdir / "variant_selection.csv",
        records,
        filtering_decisions,
        xtb_prefilter_decisions,
        stereo_reduction,
    )
    write_crest_job_plan(outdir / "crest_job_plan.csv", records, stereo_reduction)
    write_disk_and_timing_tables(outdir)


def write_enumeration_counts(path: Path, records: list[AnyLineageRecord]) -> None:
    counts: Counter[tuple[str, str, str]] = Counter(
        (record.input_molecule_id, record.molname, record.stage_name) for record in records
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["input_id", "molname", "stage", "count"],
        )
        writer.writeheader()
        for (input_id, molname, stage), count in sorted(counts.items()):
            writer.writerow(
                {"input_id": input_id, "molname": molname, "stage": stage, "count": count}
            )


def write_variant_selection(
    path: Path,
    records: list[AnyLineageRecord],
    filtering_decisions: list[FilteringDecision],
    xtb_prefilter_decisions: list[XtbPrefilterDecision],
    stereo_reduction: StereoReductionResult,
) -> None:
    record_by_id = {record.id: record for record in records}
    filtering_by_id: dict[str, list[FilteringDecision]] = defaultdict(list)
    for decision in filtering_decisions:
        filtering_by_id[decision.record_id].append(decision)
    xtb_by_stereo = {decision.stereo_id: decision for decision in xtb_prefilter_decisions}
    crest_by_stereo = {decision.stereo_id: decision for decision in stereo_reduction.decisions}
    columns = [
        "molname",
        "input_id",
        "protomer_id",
        "tautomer_id",
        "stereo_id",
        "formula",
        "charge",
        "svp_score_total",
        "protomer_penalty",
        "tautomer_penalty",
        "stereo_penalty",
        "chemistry_sanity_penalty",
        "complexity_penalty",
        "cheap_3d_energy_penalty",
        "accepted_for_3d",
        "accepted_for_xtb_prefilter",
        "accepted_for_crest",
        "accepted_for_thermo",
        "accepted_for_censo",
        "rejection_stage",
        "rejection_reason",
        "rescue_rule",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in sorted(
            [record for record in records if record.stage_name == "stereo"],
            key=lambda item: item.id,
        ):
            decisions = filtering_by_id.get(record.id, [])
            latest = decisions[-1] if decisions else None
            xtb = xtb_by_stereo.get(record.id)
            crest = crest_by_stereo.get(record.id)
            rejection_stage = ""
            rejection_reason = ""
            for stage_name, decision in (
                ("3d", latest),
                ("xtb_prefilter", xtb),
                ("crest", crest),
            ):
                if decision is not None and not getattr(decision, "selected", True):
                    rejection_stage = stage_name
                    rejection_reason = getattr(decision, "reason", "")
                    break
            writer.writerow(
                {
                    "molname": record.molname,
                    "input_id": record.input_molecule_id,
                    "protomer_id": _ancestor(record.id, "p"),
                    "tautomer_id": _ancestor(record.id, "t"),
                    "stereo_id": record.id,
                    "formula": record.molecular_formula,
                    "charge": record.formal_charge,
                    "svp_score_total": latest.score if latest else None,
                    "protomer_penalty": latest.protomer_penalty if latest else None,
                    "tautomer_penalty": latest.tautomer_penalty if latest else None,
                    "stereo_penalty": latest.stereo_penalty if latest else None,
                    "chemistry_sanity_penalty": latest.chemistry_sanity_penalty if latest else None,
                    "complexity_penalty": latest.complexity_penalty if latest else None,
                    "cheap_3d_energy_penalty": latest.cheap_3d_energy_penalty if latest else None,
                    "accepted_for_3d": latest.selected if latest else True,
                    "accepted_for_xtb_prefilter": xtb.selected if xtb else True,
                    "accepted_for_crest": crest.selected_for_crest if crest else False,
                    "accepted_for_thermo": _has_child_stage(record.id, record_by_id, "thermo"),
                    "accepted_for_censo": False,
                    "rejection_stage": rejection_stage,
                    "rejection_reason": rejection_reason,
                    "rescue_rule": latest.rescue_reason if latest else None,
                    "warnings": " | ".join(record.warnings),
                }
            )


def write_crest_job_plan(
    path: Path,
    records: list[AnyLineageRecord],
    stereo_reduction: StereoReductionResult,
) -> None:
    seed_records = {
        record.id: record for record in records if record.stage_name == "seed_conformer"
    }
    columns = [
        "molname",
        "variant_id",
        "seed_id",
        "priority_rank",
        "estimated_rotatable_bonds",
        "estimated_atoms",
        "selected_for_crest",
        "reason",
        "workdir",
        "command",
        "status",
        "elapsed_seconds",
        "disk_mb",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for rank, decision in enumerate(stereo_reduction.decisions, start=1):
            seed = seed_records.get(decision.seed_id)
            mol = getattr(seed, "rdkit_mol", None) if seed is not None else None
            writer.writerow(
                {
                    "molname": seed.molname if seed else "",
                    "variant_id": decision.stereo_id,
                    "seed_id": decision.seed_id,
                    "priority_rank": rank,
                    "estimated_rotatable_bonds": Lipinski.NumRotatableBonds(mol)
                    if mol is not None
                    else None,
                    "estimated_atoms": mol.GetNumAtoms() if mol is not None else None,
                    "selected_for_crest": decision.selected_for_crest,
                    "reason": decision.reason,
                    "workdir": "",
                    "command": "",
                    "status": "planned" if decision.selected_for_crest else "skipped",
                    "elapsed_seconds": "",
                    "disk_mb": "",
                }
            )


def write_disk_and_timing_tables(outdir: Path) -> None:
    summary = outdir / "stage_summary.csv"
    disk_path = outdir / "disk_usage_by_stage.csv"
    timing_path = outdir / "timing_by_stage.csv"
    if not summary.exists():
        disk_path.write_text("stage,run_dir_size_mb,xyz_file_count\n", encoding="utf-8")
        timing_path.write_text("stage,elapsed_seconds,status\n", encoding="utf-8")
        return
    with summary.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    with disk_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stage", "run_dir_size_mb", "xyz_file_count"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "stage": row.get("stage"),
                    "run_dir_size_mb": row.get("run_dir_size_mb"),
                    "xyz_file_count": row.get("xyz_file_count"),
                }
            )
    with timing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stage", "elapsed_seconds", "status"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "stage": row.get("stage"),
                    "elapsed_seconds": row.get("elapsed_seconds"),
                    "status": row.get("status"),
                }
            )


def _ancestor(record_id: str, marker: str) -> str | None:
    token = f"_{marker}"
    if token not in record_id:
        return None
    parts = record_id.split("_")
    for index, part in enumerate(parts):
        if part.startswith(marker) and len(part) >= 2 and part[1:3].isdigit():
            return "_".join(parts[: index + 2]) if index + 1 < len(parts) else record_id
    return None


def _has_child_stage(
    parent_fragment: str,
    record_by_id: dict[str, AnyLineageRecord],
    stage: str,
) -> bool:
    return any(
        stage == record.stage_name and parent_fragment in (record.parent_id or "")
        for record in record_by_id.values()
    )
