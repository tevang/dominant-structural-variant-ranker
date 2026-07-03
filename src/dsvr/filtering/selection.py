from __future__ import annotations

import csv
import json
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from dsvr.config import RunConfig
from dsvr.filtering.budgets import budget, seed_budget
from dsvr.filtering.variant_score import (
    PenaltyBreakdown,
    cheap_variant_score,
    score_variant,
)
from dsvr.models import SeedConformerRecord, StereoRecord


@dataclass(frozen=True)
class FilteringDecision:
    record_id: str
    input_molecule_id: str
    stage: str
    selected: bool
    rank: int | None
    score: float | None
    reason: str
    mode: str
    rescue_reason: str | None = None
    protomer_penalty: float | None = None
    tautomer_penalty: float | None = None
    stereo_penalty: float | None = None
    chemistry_sanity_penalty: float | None = None
    complexity_penalty: float | None = None
    cheap_3d_energy_penalty: float | None = None
    warnings: list[str] | None = None


def select_stereo_records(
    records: list[StereoRecord],
    config: RunConfig,
    stage: str,
) -> tuple[list[StereoRecord], list[FilteringDecision]]:
    if not records or not config.variant_filtering.enabled:
        return records, _pass_decisions(records, config, stage, "filtering_disabled")
    if config.variant_filtering.mode == "exhaustive":
        return records, _pass_decisions(records, config, stage, "exhaustive_mode_expensive")

    max_count = budget(config, "max_variants_before_3d_per_molecule")
    if stage == "cheap_score":
        max_count = budget(config, "max_variants_after_cheap_score_per_molecule")
    if max_count is None:
        return records, _pass_decisions(records, config, stage, "unlimited_budget")

    selected: list[StereoRecord] = []
    decisions: list[FilteringDecision] = []
    for input_id, group in _group_by_input(records).items():
        ranked = sorted(
            (
                (
                    score_variant(
                        record,
                        {
                            "ph": config.chemistry.ph,
                            "achiral_environment": True,
                        },
                    ),
                    record,
                )
                for record in group
            ),
            key=lambda item: (item[0].total, item[1].id),
        )
        score_by_id = {record.id: score for score, record in ranked}
        keep_ids = _rescue_reasons([record for _, record in ranked], config, score_by_id)
        ranked = sorted(
            ranked,
            key=lambda item: (
                0 if item[1].id in keep_ids else 1,
                item[0].total,
                item[1].id,
            ),
        )
        best_penalty = ranked[0][0].total if ranked else 0.0
        selected_count = 0
        for index, (breakdown, record) in enumerate(ranked, start=1):
            within_budget = selected_count < max_count
            within_cutoff = (
                breakdown.total <= config.variant_filtering.absolute_penalty_cutoff
                and breakdown.total - best_penalty
                <= config.variant_filtering.relative_penalty_cutoff
            )
            rescue_reason = keep_ids.get(record.id)
            rescued = rescue_reason is not None
            keep = within_budget and (rescued or within_cutoff)
            if keep:
                selected.append(record)
                selected_count += 1
            decisions.append(
                FilteringDecision(
                    record_id=record.id,
                    input_molecule_id=input_id,
                    stage=stage,
                    selected=keep,
                    rank=index,
                    score=breakdown.total,
                    reason=_reason(stage, keep, within_budget, within_cutoff, rescued),
                    mode=config.variant_filtering.mode,
                    rescue_reason=rescue_reason,
                    **_breakdown_fields(breakdown),
                )
            )
    return selected, decisions


def select_seed_records(
    records: list[SeedConformerRecord],
    config: RunConfig,
) -> tuple[list[SeedConformerRecord], list[FilteringDecision]]:
    if not records or not config.variant_filtering.enabled:
        return records, _pass_seed_decisions(records, config, "filtering_disabled")
    if config.variant_filtering.mode == "exhaustive":
        return records, _pass_seed_decisions(records, config, "exhaustive_mode_expensive")

    per_variant = seed_budget(config)
    per_molecule = budget(config, "max_variants_for_crest_per_molecule")
    if per_molecule is None:
        per_molecule = config.crest.max_jobs_per_molecule
    else:
        per_molecule = min(per_molecule, config.crest.max_jobs_per_molecule)
    selected: list[SeedConformerRecord] = []
    decisions: list[FilteringDecision] = []

    for input_id, group in _group_by_input(records).items():
        by_parent: dict[str | None, list[SeedConformerRecord]] = {}
        for record in group:
            by_parent.setdefault(record.parent_id, []).append(record)
        parent_ranked = sorted(
            by_parent.items(),
            key=lambda item: (_best_seed_energy(item[1]), str(item[0])),
        )
        selected_parents = {parent for parent, _ in parent_ranked[:per_molecule]}
        for parent, seeds in parent_ranked:
            ranked_seeds = sorted(
                seeds,
                key=lambda seed: (
                    float("inf") if seed.energy_kcal_mol is None else seed.energy_kcal_mol,
                    seed.id,
                ),
            )
            for index, seed in enumerate(ranked_seeds, start=1):
                keep = parent in selected_parents and (per_variant is None or index <= per_variant)
                if keep:
                    selected.append(seed)
                decisions.append(
                    FilteringDecision(
                        record_id=seed.id,
                        input_molecule_id=input_id,
                        stage="crest_seed_selection",
                        selected=keep,
                        rank=index,
                        score=seed.energy_kcal_mol,
                        reason="selected_for_crest" if keep else "over_seed_or_variant_budget",
                        mode=config.variant_filtering.mode,
                    )
                )
    return selected, decisions


def write_filtering_decisions(path: Path, decisions: list[FilteringDecision]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for decision in decisions:
            handle.write(json.dumps(asdict(decision), sort_keys=True) + "\n")


def write_filtering_csv(path: Path, decisions: list[FilteringDecision]) -> None:
    columns = [
        "record_id",
        "input_molecule_id",
        "stage",
        "selected",
        "rank",
        "score",
        "reason",
        "mode",
        "rescue_reason",
        "protomer_penalty",
        "tautomer_penalty",
        "stereo_penalty",
        "chemistry_sanity_penalty",
        "complexity_penalty",
        "cheap_3d_energy_penalty",
        "warnings",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for decision in decisions:
            row = asdict(decision)
            if isinstance(row.get("warnings"), list):
                row["warnings"] = " | ".join(row["warnings"])
            writer.writerow(row)


def write_penalty_outputs(path: Path, decisions: list[FilteringDecision]) -> None:
    path.mkdir(parents=True, exist_ok=True)
    penalty_decisions = [decision for decision in decisions if decision.score is not None]
    _write_decision_subset(path / "variant_penalties.csv", penalty_decisions)
    _write_decision_subset(
        path / "accepted_variants.csv",
        [decision for decision in penalty_decisions if decision.selected],
    )
    _write_decision_subset(
        path / "rejected_variants.csv",
        [decision for decision in penalty_decisions if not decision.selected],
    )
    with (path / "penalty_breakdown.jsonl").open("w", encoding="utf-8") as handle:
        for decision in penalty_decisions:
            handle.write(json.dumps(asdict(decision), sort_keys=True) + "\n")


def _group_by_input(records: Iterable[StereoRecord | SeedConformerRecord]):
    grouped: dict[str, list[StereoRecord | SeedConformerRecord]] = {}
    for record in records:
        grouped.setdefault(record.input_molecule_id, []).append(record)
    return grouped


def _rescue_reasons(
    records: list[StereoRecord],
    config: RunConfig,
    score_by_id: dict[str, PenaltyBreakdown],
) -> dict[str, str]:
    if not config.variant_filtering.rescue_rules_enabled:
        return {}
    rescued: dict[str, str] = {}
    if config.variant_filtering.keep_original_state and records:
        rescued[records[0].id] = "rescue_original_input_state"
    if config.variant_filtering.keep_best_per_charge_state:
        _add_rescue(
            rescued,
            _best_per_key(records, lambda record: record.formal_charge, score_by_id),
            "rescue_best_per_formal_charge",
        )
    if config.variant_filtering.keep_best_per_formula:
        _add_rescue(
            rescued,
            _best_per_key(records, lambda record: record.molecular_formula, score_by_id),
            "rescue_best_per_formula_or_proton_count",
        )
    if config.variant_filtering.keep_best_per_protomer:
        _add_rescue(
            rescued,
            _best_per_key(records, lambda record: _ancestor(record.id, "_p"), score_by_id),
            "rescue_best_per_parent_protomer",
        )
    if config.variant_filtering.keep_best_per_tautomer_family:
        _add_rescue(
            rescued,
            _best_per_key(records, lambda record: _ancestor(record.id, "_t"), score_by_id),
            "rescue_best_per_tautomer_family",
        )
    return rescued


def _add_rescue(rescued: dict[str, str], ids: set[str], reason: str) -> None:
    for record_id in ids:
        rescued.setdefault(record_id, reason)


def _best_per_key(
    records: list[StereoRecord],
    key_fn,
    score_by_id: dict[str, PenaltyBreakdown],
) -> set[str]:
    best: dict[object, tuple[float, str]] = {}
    for record in records:
        key = key_fn(record)
        score = score_by_id.get(record.id)
        value = score.total if score is not None else cheap_variant_score(record).penalty
        current = best.get(key)
        if current is None or (value, record.id) < current:
            best[key] = (value, record.id)
    return {record_id for _, record_id in best.values()}


def _ancestor(record_id: str, marker: str) -> str:
    index = record_id.find(marker)
    if index < 0:
        return record_id
    next_marker = record_id.find("_", index + len(marker) + 2)
    return record_id if next_marker < 0 else record_id[:next_marker]


def _best_seed_energy(records: list[SeedConformerRecord]) -> float:
    values = [record.energy_kcal_mol for record in records if record.energy_kcal_mol is not None]
    return min(values) if values else float("inf")


def _reason(
    stage: str,
    selected: bool,
    within_budget: bool,
    within_cutoff: bool,
    rescued: bool,
) -> str:
    if rescued:
        return f"{stage}_rescue_rule"
    if selected:
        return f"{stage}_selected"
    if not within_budget:
        return f"{stage}_over_budget"
    if not within_cutoff:
        return f"{stage}_over_penalty_cutoff"
    return f"{stage}_filtered"


def _breakdown_fields(breakdown: PenaltyBreakdown) -> dict[str, object]:
    return {
        "protomer_penalty": breakdown.protomer_penalty,
        "tautomer_penalty": breakdown.tautomer_penalty,
        "stereo_penalty": breakdown.stereo_penalty,
        "chemistry_sanity_penalty": breakdown.chemistry_sanity_penalty,
        "complexity_penalty": breakdown.complexity_penalty,
        "cheap_3d_energy_penalty": breakdown.cheap_3d_energy_penalty,
        "warnings": breakdown.warnings,
    }


def _write_decision_subset(path: Path, decisions: list[FilteringDecision]) -> None:
    columns = [
        "record_id",
        "input_molecule_id",
        "stage",
        "selected",
        "rank",
        "score",
        "reason",
        "mode",
        "rescue_reason",
        "protomer_penalty",
        "tautomer_penalty",
        "stereo_penalty",
        "chemistry_sanity_penalty",
        "complexity_penalty",
        "cheap_3d_energy_penalty",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for decision in decisions:
            row = asdict(decision)
            if isinstance(row.get("warnings"), list):
                row["warnings"] = " | ".join(row["warnings"])
            writer.writerow({column: row.get(column) for column in columns})


def _pass_decisions(
    records: list[StereoRecord],
    config: RunConfig,
    stage: str,
    reason: str,
) -> list[FilteringDecision]:
    return [
        FilteringDecision(
            record_id=record.id,
            input_molecule_id=record.input_molecule_id,
            stage=stage,
            selected=True,
            rank=index,
            score=None,
            reason=reason,
            mode=config.variant_filtering.mode,
        )
        for index, record in enumerate(records, start=1)
    ]


def _pass_seed_decisions(
    records: list[SeedConformerRecord],
    config: RunConfig,
    reason: str,
) -> list[FilteringDecision]:
    return [
        FilteringDecision(
            record_id=record.id,
            input_molecule_id=record.input_molecule_id,
            stage="crest_seed_selection",
            selected=True,
            rank=index,
            score=record.energy_kcal_mol,
            reason=reason,
            mode=config.variant_filtering.mode,
        )
        for index, record in enumerate(records, start=1)
    ]
