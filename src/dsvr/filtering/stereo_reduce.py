from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.models import (
    CrestConformerRecord,
    SeedConformerRecord,
    StereoRecord,
    ThermoRecord,
    short_hash,
)


@dataclass(frozen=True)
class StereoReductionDecision:
    seed_id: str
    stereo_id: str | None
    input_molecule_id: str
    selected_for_crest: bool
    representative_seed_id: str | None
    representative_stereo_id: str | None
    relationship: str
    reason: str
    group_key: str


@dataclass(frozen=True)
class StereoReductionResult:
    selected_seeds: list[SeedConformerRecord]
    decisions: list[StereoReductionDecision]
    representative_to_equivalent_seed_ids: dict[str, list[str]]
    jobs_saved: int


def reduce_seeds_for_crest(
    seeds: list[SeedConformerRecord],
    stereo_records: list[StereoRecord],
    config: RunConfig,
) -> StereoReductionResult:
    if not _enabled(config) or not seeds:
        disabled_decisions = [
            StereoReductionDecision(
                seed_id=seed.id,
                stereo_id=seed.parent_id,
                input_molecule_id=seed.input_molecule_id,
                selected_for_crest=True,
                representative_seed_id=seed.id,
                representative_stereo_id=seed.parent_id,
                relationship="not_reduced",
                reason="stereo_reduction_disabled",
                group_key=_seed_group_key(seed),
            )
            for seed in seeds
        ]
        return StereoReductionResult(seeds, disabled_decisions, {}, jobs_saved=0)

    stereo_by_id = {record.id: record for record in stereo_records}
    enantiomer_representatives = _enantiomer_representatives(stereo_records)
    selected: list[SeedConformerRecord] = []
    decisions: list[StereoReductionDecision] = []
    representative_to_equivalent_seed_ids: dict[str, list[str]] = defaultdict(list)

    by_stereo: dict[str | None, list[SeedConformerRecord]] = defaultdict(list)
    for seed in seeds:
        by_stereo[seed.parent_id].append(seed)

    representative_stereo_ids = set(enantiomer_representatives.values())
    skipped_stereo_ids = set(enantiomer_representatives) - representative_stereo_ids
    for stereo_id, stereo_seeds in sorted(by_stereo.items(), key=lambda item: str(item[0])):
        representative_stereo_id = enantiomer_representatives.get(stereo_id or "")
        relationship = "enantiomer_pair" if representative_stereo_id else "not_reduced"
        if representative_stereo_id is None or representative_stereo_id == stereo_id:
            selected.extend(stereo_seeds)
            for seed in stereo_seeds:
                decisions.append(
                    StereoReductionDecision(
                        seed_id=seed.id,
                        stereo_id=stereo_id,
                        input_molecule_id=seed.input_molecule_id,
                        selected_for_crest=True,
                        representative_seed_id=seed.id,
                        representative_stereo_id=stereo_id,
                        relationship=relationship,
                        reason=(
                            "representative_enantiomer_seed_selected_for_crest"
                            if stereo_id in representative_stereo_ids
                            else "selected_for_crest"
                        ),
                        group_key=_stereo_group_key(stereo_by_id.get(stereo_id or ""), seed),
                    )
                )
            continue

        representative_seeds = by_stereo.get(representative_stereo_id, [])
        for seed in stereo_seeds:
            representative_seed = _matching_representative_seed(seed, representative_seeds)
            if representative_seed is not None:
                representative_to_equivalent_seed_ids[representative_seed.id].append(seed.id)
            decisions.append(
                StereoReductionDecision(
                    seed_id=seed.id,
                    stereo_id=stereo_id,
                    input_molecule_id=seed.input_molecule_id,
                    selected_for_crest=False,
                    representative_seed_id=representative_seed.id if representative_seed else None,
                    representative_stereo_id=representative_stereo_id,
                    relationship="enantiomer_pair",
                    reason=(
                        "crest_skipped_for_enantiomer_equivalent_in_achiral_solvent"
                        if stereo_id in skipped_stereo_ids
                        else "crest_skipped_by_stereo_reduction"
                    ),
                    group_key=_stereo_group_key(stereo_by_id.get(stereo_id or ""), seed),
                )
            )

    jobs_saved = sum(1 for decision in decisions if not decision.selected_for_crest)
    return StereoReductionResult(
        selected_seeds=selected,
        decisions=decisions,
        representative_to_equivalent_seed_ids=dict(representative_to_equivalent_seed_ids),
        jobs_saved=jobs_saved,
    )


def expand_enantiomer_mapped_crest_records(
    records: list[CrestConformerRecord],
    seed_by_id: dict[str, SeedConformerRecord],
    reduction: StereoReductionResult,
    config: RunConfig,
) -> list[CrestConformerRecord]:
    if not config.stereo_filtering.keep_mapping_to_all_stereo_outputs:
        return records
    expanded = list(records)
    for record in records:
        representative_seed_id = record.parent_id
        if representative_seed_id is None:
            continue
        for equivalent_seed_id in reduction.representative_to_equivalent_seed_ids.get(
            representative_seed_id,
            [],
        ):
            equivalent_seed = seed_by_id.get(equivalent_seed_id)
            if equivalent_seed is None:
                continue
            expanded.append(_copy_crest_record_for_equivalent_seed(record, equivalent_seed))
    return expanded


def expand_enantiomer_mapped_thermo_records(
    records: list[ThermoRecord],
    expanded_crest_records: list[CrestConformerRecord],
    config: RunConfig,
) -> list[ThermoRecord]:
    if not config.stereo_filtering.keep_mapping_to_all_stereo_outputs:
        return records
    equivalent_crest_by_representative: dict[str, list[CrestConformerRecord]] = defaultdict(list)
    for crest_record in expanded_crest_records:
        stereo_reduction = crest_record.metadata.get("stereo_reduction")
        if not isinstance(stereo_reduction, dict):
            continue
        representative_id = stereo_reduction.get("mapped_from_representative_conformer_id")
        if isinstance(representative_id, str):
            equivalent_crest_by_representative[representative_id].append(crest_record)

    expanded = list(records)
    for record in records:
        representative_crest_id = record.parent_id
        if representative_crest_id is None:
            continue
        for equivalent_crest in equivalent_crest_by_representative.get(representative_crest_id, []):
            expanded.append(_copy_thermo_record_for_equivalent_crest(record, equivalent_crest))
    return expanded


def write_stereo_reduction_outputs(path: Path, reduction: StereoReductionResult) -> None:
    path.mkdir(parents=True, exist_ok=True)
    rows = [asdict(decision) for decision in reduction.decisions]
    columns = [
        "seed_id",
        "stereo_id",
        "input_molecule_id",
        "selected_for_crest",
        "representative_seed_id",
        "representative_stereo_id",
        "relationship",
        "reason",
        "group_key",
    ]
    with (path / "stereo_reduction.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    with (path / "stereo_reduction.jsonl").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def _enabled(config: RunConfig) -> bool:
    return (
        config.stereo_filtering.collapse_enantiomers_in_achiral_solvent
        and not config.stereo_filtering.solvent_is_chiral
        and config.stereo_filtering.run_crest_for_enantiomer_pairs_once
    )


def _enantiomer_representatives(records: list[StereoRecord]) -> dict[str, str]:
    groups: dict[tuple[str | None, int | None, str], list[StereoRecord]] = defaultdict(list)
    for record in records:
        groups[_stereo_group_tuple(record)].append(record)

    representatives: dict[str, str] = {}
    for group in groups.values():
        if len(group) != 2 or not all(_single_assigned_chiral_center(record) for record in group):
            continue
        if not _opposite_single_center_configuration(group[0], group[1]):
            continue
        representative = sorted(group, key=lambda item: item.id)[0]
        for record in group:
            representatives[record.id] = representative.id
    return representatives


def _stereo_group_key(
    record: StereoRecord | None,
    fallback_seed: SeedConformerRecord | None = None,
) -> str:
    source = record or fallback_seed
    if source is None:
        return ""
    key = _source_group_tuple(source)
    return "|".join("" if value is None else str(value) for value in key)


def _stereo_group_tuple(record: StereoRecord) -> tuple[str | None, int | None, str]:
    return _source_group_tuple(record)


def _source_group_tuple(
    source: StereoRecord | SeedConformerRecord,
) -> tuple[str | None, int | None, str]:
    mol = _mol(source)
    if mol is None:
        achiral_smiles = source.canonical_smiles or source.isomeric_smiles or ""
    else:
        achiral_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    return (source.molecular_formula, source.formal_charge, achiral_smiles)


def _seed_group_key(seed: SeedConformerRecord) -> str:
    key = _stereo_group_key(None, seed)
    return str(key or "")


def _mol(record: StereoRecord | SeedConformerRecord) -> Chem.Mol | None:
    if record.rdkit_mol is not None:
        return record.rdkit_mol
    if record.isomeric_smiles:
        return Chem.MolFromSmiles(record.isomeric_smiles)
    if record.canonical_smiles:
        return Chem.MolFromSmiles(record.canonical_smiles)
    return None


def _single_assigned_chiral_center(record: StereoRecord) -> bool:
    mol = _mol(record)
    if mol is None:
        return False
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    return len(centers) == 1


def _opposite_single_center_configuration(first: StereoRecord, second: StereoRecord) -> bool:
    first_centers = _assigned_chiral_centers(first)
    second_centers = _assigned_chiral_centers(second)
    if len(first_centers) != 1 or len(second_centers) != 1:
        return False
    first_index, first_label = first_centers[0]
    second_index, second_label = second_centers[0]
    return first_index == second_index and {first_label, second_label} == {"R", "S"}


def _assigned_chiral_centers(record: StereoRecord) -> list[tuple[int, str]]:
    mol = _mol(record)
    if mol is None:
        return []
    return Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)


def _matching_representative_seed(
    seed: SeedConformerRecord,
    representative_seeds: list[SeedConformerRecord],
) -> SeedConformerRecord | None:
    if not representative_seeds:
        return None
    for representative_seed in representative_seeds:
        if representative_seed.conformer_index == seed.conformer_index:
            return representative_seed
    return sorted(
        representative_seeds,
        key=lambda item: (
            float("inf") if item.energy_kcal_mol is None else item.energy_kcal_mol,
            item.id,
        ),
    )[0]


def _copy_crest_record_for_equivalent_seed(
    record: CrestConformerRecord,
    equivalent_seed: SeedConformerRecord,
) -> CrestConformerRecord:
    metadata = dict(record.metadata)
    metadata["stereo_reduction"] = {
        "mapped_from_representative_conformer_id": record.id,
        "representative_seed_id": record.parent_id,
        "equivalent_seed_id": equivalent_seed.id,
        "assumption": "enantiomers have identical energy in achiral solvent",
    }
    new_id = f"{record.id}_stereo_equiv_{short_hash(equivalent_seed.id)}"
    return record.model_copy(
        update={
            "id": new_id,
            "parent_id": equivalent_seed.id,
            "canonical_smiles": equivalent_seed.canonical_smiles,
            "isomeric_smiles": equivalent_seed.isomeric_smiles,
            "molecular_formula": equivalent_seed.molecular_formula,
            "formal_charge": equivalent_seed.formal_charge,
            "explicit_proton_count": equivalent_seed.explicit_proton_count,
            "warnings": sorted(
                {
                    *record.warnings,
                    (
                        "Energy/population source was mapped from an enantiomeric "
                        "representative because the solvent is treated as achiral."
                    ),
                }
            ),
            "metadata": metadata,
        }
    )


def _copy_thermo_record_for_equivalent_crest(
    record: ThermoRecord,
    equivalent_crest: CrestConformerRecord,
) -> ThermoRecord:
    metadata = dict(record.metadata)
    metadata["stereo_reduction"] = {
        "mapped_from_representative_thermo_id": record.id,
        "representative_conformer_id": record.parent_id,
        "equivalent_conformer_id": equivalent_crest.id,
        "assumption": "enantiomers have identical energy in achiral solvent",
    }
    new_id = f"{record.id}_stereo_equiv_{short_hash(equivalent_crest.id)}"
    return record.model_copy(
        update={
            "id": new_id,
            "parent_id": equivalent_crest.id,
            "canonical_smiles": equivalent_crest.canonical_smiles,
            "isomeric_smiles": equivalent_crest.isomeric_smiles,
            "molecular_formula": equivalent_crest.molecular_formula,
            "formal_charge": equivalent_crest.formal_charge,
            "explicit_proton_count": equivalent_crest.explicit_proton_count,
            "warnings": sorted(
                {
                    *record.warnings,
                    (
                        "Thermo values were mapped from an enantiomeric representative "
                        "because the solvent is treated as achiral."
                    ),
                }
            ),
            "metadata": metadata,
        }
    )
