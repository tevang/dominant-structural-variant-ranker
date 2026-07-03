from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

from dsvr.config import RunConfig
from dsvr.models import SeedConformerRecord
from dsvr.runners.xtb_prefilter_runner import XtbPrefilterResult, run_xtb_prefilter


@dataclass(frozen=True)
class XtbPrefilterDecision:
    seed_id: str
    stereo_id: str | None
    input_molecule_id: str
    molname: str
    formula: str | None
    charge: int | None
    proton_count: int | None
    energy_kcal_mol: float | None
    relative_energy_kcal_mol: float | None
    selected: bool
    reason: str
    workdir: str | None = None
    warnings: list[str] | None = None


def apply_xtb_prefilter(
    seeds: list[SeedConformerRecord],
    config: RunConfig,
) -> tuple[list[SeedConformerRecord], list[XtbPrefilterDecision]]:
    if not config.xtb_prefilter.enabled or not seeds:
        disabled_decisions = [
            XtbPrefilterDecision(
                seed_id=seed.id,
                stereo_id=seed.parent_id,
                input_molecule_id=seed.input_molecule_id,
                molname=seed.molname,
                formula=seed.molecular_formula,
                charge=seed.formal_charge,
                proton_count=seed.explicit_proton_count,
                energy_kcal_mol=seed.energy_kcal_mol,
                relative_energy_kcal_mol=None,
                selected=True,
                reason="xtb_prefilter_disabled",
            )
            for seed in seeds
        ]
        return seeds, disabled_decisions

    best_seed_by_variant = _best_seed_by_variant(seeds)
    results = {
        seed.id: run_xtb_prefilter(seed, config)
        for seed in best_seed_by_variant.values()
    }
    selected_variant_ids, relative_by_seed = _selected_variant_ids(
        best_seed_by_variant,
        results,
        config,
    )
    selected_seeds = [seed for seed in seeds if seed.parent_id in selected_variant_ids]
    decisions: list[XtbPrefilterDecision] = []
    for seed in seeds:
        representative = best_seed_by_variant.get(seed.parent_id)
        result = results.get(representative.id) if representative is not None else None
        selected = seed.parent_id in selected_variant_ids
        decisions.append(
            XtbPrefilterDecision(
                seed_id=seed.id,
                stereo_id=seed.parent_id,
                input_molecule_id=seed.input_molecule_id,
                molname=seed.molname,
                formula=seed.molecular_formula,
                charge=seed.formal_charge,
                proton_count=seed.explicit_proton_count,
                energy_kcal_mol=result.energy_kcal_mol if result else None,
                relative_energy_kcal_mol=relative_by_seed.get(representative.id)
                if representative
                else None,
                selected=selected,
                reason="selected_xtb_prefilter_survivor"
                if selected
                else "rejected_by_xtb_prefilter_energy_or_budget",
                workdir=str(result.workdir) if result else None,
                warnings=result.warnings if result else [],
            )
        )
    return selected_seeds, decisions


def write_xtb_prefilter_outputs(path: Path, decisions: list[XtbPrefilterDecision]) -> None:
    path.mkdir(parents=True, exist_ok=True)
    columns = [
        "seed_id",
        "stereo_id",
        "input_molecule_id",
        "molname",
        "formula",
        "charge",
        "proton_count",
        "energy_kcal_mol",
        "relative_energy_kcal_mol",
        "selected",
        "reason",
        "workdir",
        "warnings",
    ]
    with (path / "xtb_prefilter_decisions.csv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for decision in decisions:
            row = asdict(decision)
            if isinstance(row.get("warnings"), list):
                row["warnings"] = " | ".join(row["warnings"])
            writer.writerow(row)
    with (path / "xtb_prefilter_decisions.jsonl").open("w", encoding="utf-8") as handle:
        for decision in decisions:
            handle.write(json.dumps(asdict(decision), sort_keys=True, default=str) + "\n")


def _best_seed_by_variant(
    seeds: list[SeedConformerRecord],
) -> dict[str | None, SeedConformerRecord]:
    best: dict[str | None, SeedConformerRecord] = {}
    for seed in seeds:
        current = best.get(seed.parent_id)
        if current is None or _seed_key(seed) < _seed_key(current):
            best[seed.parent_id] = seed
    return best


def _seed_key(seed: SeedConformerRecord) -> tuple[float, str]:
    return (
        float("inf") if seed.energy_kcal_mol is None else seed.energy_kcal_mol,
        seed.id,
    )


def _selected_variant_ids(
    best_seed_by_variant: dict[str | None, SeedConformerRecord],
    results: dict[str, XtbPrefilterResult],
    config: RunConfig,
) -> tuple[set[str | None], dict[str, float]]:
    by_input: dict[str, list[SeedConformerRecord]] = defaultdict(list)
    for seed in best_seed_by_variant.values():
        by_input[seed.input_molecule_id].append(seed)

    selected: set[str | None] = set()
    relative_by_seed: dict[str, float] = {}
    for group in by_input.values():
        finite = [
            seed
            for seed in group
            if results.get(seed.id) is not None and results[seed.id].energy_kcal_mol is not None
        ]
        if not finite:
            selected.update(
                seed.parent_id
                for seed in group[: config.xtb_prefilter.keep_top_n_per_molecule]
            )
            continue
        best_energy = min(
            energy
            for seed in finite
            for energy in [results[seed.id].energy_kcal_mol]
            if energy is not None
        )
        ranked = sorted(
            finite,
            key=lambda seed: (results[seed.id].energy_kcal_mol, seed.id),
        )
        for seed in ranked:
            energy = results[seed.id].energy_kcal_mol
            relative = 0.0 if energy is None else energy - (best_energy or energy)
            relative_by_seed[seed.id] = relative
            if (
                relative <= config.xtb_prefilter.keep_within_kcal_mol
                and len([sid for sid in selected if _same_input(sid, best_seed_by_variant, seed)])
                < config.xtb_prefilter.keep_top_n_per_molecule
            ):
                selected.add(seed.parent_id)
        selected.update(
            seed.parent_id
            for seed in ranked[: config.xtb_prefilter.keep_top_n_per_molecule]
        )
        selected.update(
            _best_per_key(
                ranked,
                results,
                lambda seed: seed.formal_charge,
                config.xtb_prefilter.keep_top_n_per_charge,
            )
        )
        selected.update(
            _best_per_key(
                ranked,
                results,
                lambda seed: (seed.molecular_formula, seed.explicit_proton_count),
                config.xtb_prefilter.keep_top_n_per_formula,
            )
        )
        selected_for_group = [seed for seed in ranked if seed.parent_id in selected]
        if len(selected_for_group) > config.xtb_prefilter.max_variants_per_molecule:
            keep = {
                seed.parent_id
                for seed in selected_for_group[: config.xtb_prefilter.max_variants_per_molecule]
            }
            selected.difference_update(seed.parent_id for seed in selected_for_group)
            selected.update(keep)
    return selected, relative_by_seed


def _same_input(
    stereo_id: str | None,
    best_seed_by_variant: dict[str | None, SeedConformerRecord],
    seed: SeedConformerRecord,
) -> bool:
    other = best_seed_by_variant.get(stereo_id)
    return other is not None and other.input_molecule_id == seed.input_molecule_id


def _best_per_key(
    seeds: list[SeedConformerRecord],
    results: dict[str, XtbPrefilterResult],
    key_fn,
    limit: int,
) -> set[str | None]:
    grouped: dict[object, list[SeedConformerRecord]] = defaultdict(list)
    for seed in seeds:
        grouped[key_fn(seed)].append(seed)
    selected: set[str | None] = set()
    for group in grouped.values():
        ranked = sorted(group, key=lambda seed: (results[seed.id].energy_kcal_mol, seed.id))
        selected.update(seed.parent_id for seed in ranked[:limit])
    return selected
