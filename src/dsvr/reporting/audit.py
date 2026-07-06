from __future__ import annotations

import csv
import shutil
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
    _publish_plausible_state_tables(outdir)
    write_variant_decisions(
        outdir / "variant_decisions.csv",
        outdir=outdir,
        records=records,
        filtering_decisions=filtering_decisions,
        xtb_prefilter_decisions=xtb_prefilter_decisions,
        stereo_reduction=stereo_reduction,
    )
    write_variant_selection(
        outdir / "variant_selection.csv",
        records,
        filtering_decisions,
        xtb_prefilter_decisions,
        stereo_reduction,
    )
    write_crest_job_plan(outdir / "crest_job_plan.csv", records, stereo_reduction)
    write_disk_and_timing_tables(outdir)


VARIANT_DECISION_COLUMNS = [
    "input_id",
    "molname",
    "protomer_id",
    "tautomer_id",
    "stereoisomer_id",
    "final_variant_id",
    "stage",
    "selected",
    "rejection_reason",
    "rescue_rule",
    "formal_charge",
    "formula",
    "tautomer_relative_energy_kcal_mol",
    "stereo_relative_energy_kcal_mol",
    "final_auto3d_energy",
    "warnings",
]


def write_variant_decisions(
    path: Path,
    *,
    outdir: Path,
    records: list[AnyLineageRecord],
    filtering_decisions: list[FilteringDecision],
    xtb_prefilter_decisions: list[XtbPrefilterDecision],
    stereo_reduction: StereoReductionResult,
) -> None:
    rows: list[dict[str, object]] = []
    rows.extend(_protomer_decision_rows(outdir))
    rows.extend(_tautomer_decision_rows(outdir))
    rows.extend(_stereo_decision_rows(outdir))
    rows.extend(
        _pipeline_decision_rows(
            records,
            filtering_decisions,
            xtb_prefilter_decisions,
            stereo_reduction,
        )
    )
    rows.extend(_final_variant_decision_rows(records))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=VARIANT_DECISION_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in VARIANT_DECISION_COLUMNS})


def _publish_plausible_state_tables(outdir: Path) -> None:
    sources = {
        "protomers_all.csv": outdir / "enumeration" / "protomers" / "protomers_all.csv",
        "protomers_selected.csv": outdir / "enumeration" / "protomers" / "protomers_selected.csv",
        "protomers_rejected.csv": outdir / "enumeration" / "protomers" / "protomers_rejected.csv",
        "tautomers_all_pre_auto3d.csv": outdir / "enumeration" / "tautomers" / "tautomers_all_pre_auto3d.csv",
        "tautomers_auto3d_ranked.csv": outdir / "enumeration" / "tautomers" / "tautomers_auto3d_ranked.csv",
        "tautomers_selected.csv": outdir / "enumeration" / "tautomers" / "tautomers_selected.csv",
        "tautomers_rejected.csv": outdir / "enumeration" / "tautomers" / "tautomers_rejected.csv",
    }
    for name, source in sources.items():
        target = outdir / name
        if source.exists() and source.resolve() != target.resolve():
            shutil.copy2(source, target)
        elif not target.exists():
            target.write_text("", encoding="utf-8")


def _protomer_decision_rows(outdir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in _read_csv(outdir / "protomers_selected.csv"):
        rows.append(
            _decision_row(
                input_id=row.get("input_molecule_id"),
                molname=row.get("molname"),
                protomer_id=row.get("protomer_id"),
                stage="protomer",
                selected=True,
                rescue_rule=_rescue_from_reason(row.get("selection_reason")),
                formal_charge=row.get("formal_charge"),
                formula=row.get("molecular_formula"),
                warnings=row.get("warnings"),
            )
        )
    for row in _read_csv(outdir / "protomers_rejected.csv"):
        rows.append(
            _decision_row(
                input_id=row.get("input_molecule_id"),
                molname=row.get("molname"),
                protomer_id=row.get("protomer_id") or row.get("candidate_index"),
                stage="protomer",
                selected=False,
                rejection_reason=row.get("selection_reason") or row.get("reason"),
                formal_charge=row.get("formal_charge"),
                formula=row.get("molecular_formula"),
                warnings=row.get("warnings"),
            )
        )
    return rows


def _tautomer_decision_rows(outdir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in _read_csv(outdir / "tautomers_auto3d_ranked.csv"):
        selected = _bool(row.get("selected"))
        rows.append(
            _decision_row(
                input_id=row.get("input_molecule_id"),
                molname=row.get("molname"),
                protomer_id=row.get("protomer_id"),
                tautomer_id=row.get("tautomer_id"),
                stage="tautomer",
                selected=selected,
                rejection_reason="" if selected else row.get("reason"),
                rescue_rule=_tautomer_rescue_rule(row),
                tautomer_relative_energy_kcal_mol=row.get("relative_energy_kcal_mol"),
                warnings=row.get("warnings"),
            )
        )
    return rows


def _stereo_decision_rows(outdir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in _read_csv(outdir / "stereoisomers_all.csv"):
        selected = _bool(row.get("selected"))
        rows.append(
            _decision_row(
                input_id=row.get("input_molecule_id"),
                molname=row.get("molname"),
                protomer_id=_ancestor(row.get("id"), "p"),
                tautomer_id=row.get("parent_id") or _ancestor(row.get("id"), "t"),
                stereoisomer_id=row.get("id"),
                stage="stereoisomer",
                selected=selected,
                rejection_reason="" if selected else row.get("reason"),
                rescue_rule="enantiomer_mapped" if row.get("relationship") == "enantiomer_mapped" else "",
                formal_charge=row.get("formal_charge"),
                formula=row.get("molecular_formula"),
                stereo_relative_energy_kcal_mol=row.get("relative_energy_kcal_mol"),
                warnings=row.get("warnings"),
            )
        )
    return rows


def _pipeline_decision_rows(
    records: list[AnyLineageRecord],
    filtering_decisions: list[FilteringDecision],
    xtb_prefilter_decisions: list[XtbPrefilterDecision],
    stereo_reduction: StereoReductionResult,
) -> list[dict[str, object]]:
    records_by_id = {record.id: record for record in records}
    rows: list[dict[str, object]] = []
    for decision in filtering_decisions:
        record = records_by_id.get(decision.record_id)
        rows.append(
            _decision_row(
                input_id=decision.input_molecule_id,
                molname=record.molname if record else "",
                protomer_id=_ancestor(decision.record_id, "p"),
                tautomer_id=_ancestor(decision.record_id, "t"),
                stereoisomer_id=_ancestor(decision.record_id, "s"),
                final_variant_id=decision.record_id if decision.stage == "crest_seed_selection" else "",
                stage=decision.stage,
                selected=decision.selected,
                rejection_reason="" if decision.selected else decision.reason,
                rescue_rule=decision.rescue_reason,
                formal_charge=record.formal_charge if record else "",
                formula=record.molecular_formula if record else "",
                final_auto3d_energy=getattr(record, "energy_kcal_mol", ""),
                warnings=" | ".join(decision.warnings or []) or _record_warnings(record),
            )
        )
    for decision in xtb_prefilter_decisions:
        rows.append(
            _decision_row(
                input_id=decision.input_molecule_id,
                molname=decision.molname,
                protomer_id=_ancestor(decision.stereo_id, "p"),
                tautomer_id=_ancestor(decision.stereo_id, "t"),
                stereoisomer_id=decision.stereo_id,
                final_variant_id=decision.seed_id,
                stage="xtb_prefilter",
                selected=decision.selected,
                rejection_reason="" if decision.selected else decision.reason,
                formal_charge=decision.charge,
                formula=decision.formula,
                final_auto3d_energy=decision.energy_kcal_mol,
                warnings=" | ".join(decision.warnings or []),
            )
        )
    for decision in stereo_reduction.decisions:
        rows.append(
            _decision_row(
                input_id=decision.input_molecule_id,
                protomer_id=_ancestor(decision.stereo_id, "p"),
                tautomer_id=_ancestor(decision.stereo_id, "t"),
                stereoisomer_id=decision.stereo_id,
                final_variant_id=decision.seed_id,
                stage="crest_stereo_reduction",
                selected=decision.selected_for_crest,
                rejection_reason="" if decision.selected_for_crest else decision.reason,
                rescue_rule="representative_seed" if decision.relationship != "not_reduced" else "",
            )
        )
    return rows


def _final_variant_decision_rows(records: list[AnyLineageRecord]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in records:
        if record.stage_name != "ranked_variant":
            continue
        rows.append(
            _decision_row(
                input_id=record.input_molecule_id,
                molname=record.molname,
                protomer_id=_ancestor(record.parent_id, "p"),
                tautomer_id=_ancestor(record.parent_id, "t"),
                stereoisomer_id=_ancestor(record.parent_id, "s"),
                final_variant_id=record.id,
                stage="final_variant",
                selected=True,
                formal_charge=record.formal_charge,
                formula=record.molecular_formula,
                final_auto3d_energy=record.score_kcal_mol if hasattr(record, "score_kcal_mol") else "",
                warnings=_record_warnings(record),
            )
        )
    return rows


def _decision_row(
    *,
    input_id: object = "",
    molname: object = "",
    protomer_id: object = "",
    tautomer_id: object = "",
    stereoisomer_id: object = "",
    final_variant_id: object = "",
    stage: object = "",
    selected: object = "",
    rejection_reason: object = "",
    rescue_rule: object = "",
    formal_charge: object = "",
    formula: object = "",
    tautomer_relative_energy_kcal_mol: object = "",
    stereo_relative_energy_kcal_mol: object = "",
    final_auto3d_energy: object = "",
    warnings: object = "",
) -> dict[str, object]:
    return {
        "input_id": input_id or "",
        "molname": molname or "",
        "protomer_id": protomer_id or "",
        "tautomer_id": tautomer_id or "",
        "stereoisomer_id": stereoisomer_id or "",
        "final_variant_id": final_variant_id or "",
        "stage": stage or "",
        "selected": selected,
        "rejection_reason": rejection_reason or "",
        "rescue_rule": rescue_rule or "",
        "formal_charge": formal_charge if formal_charge is not None else "",
        "formula": formula or "",
        "tautomer_relative_energy_kcal_mol": tautomer_relative_energy_kcal_mol
        if tautomer_relative_energy_kcal_mol is not None
        else "",
        "stereo_relative_energy_kcal_mol": stereo_relative_energy_kcal_mol
        if stereo_relative_energy_kcal_mol is not None
        else "",
        "final_auto3d_energy": final_auto3d_energy if final_auto3d_energy is not None else "",
        "warnings": warnings or "",
    }


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def _rescue_from_reason(reason: object) -> str:
    text = str(reason or "")
    if text in {"input_reference_state", "best_representative_per_charge"}:
        return text
    return ""


def _tautomer_rescue_rule(row: dict[str, str]) -> str:
    reason = row.get("reason", "")
    warnings = row.get("warnings", "")
    if _bool(row.get("selected")) and "fallback" in reason.lower():
        return reason
    if "TAUTOMER_TIMEOUT_FALLBACK" in warnings:
        return "TAUTOMER_TIMEOUT_FALLBACK"
    return ""


def _record_warnings(record: AnyLineageRecord | None) -> str:
    return " | ".join(getattr(record, "warnings", []) or [])


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


def _ancestor(record_id: str | None, marker: str) -> str | None:
    if record_id is None:
        return None
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
