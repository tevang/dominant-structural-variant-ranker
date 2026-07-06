from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.models import StereoRecord
from dsvr.runners.auto3d_runner import Auto3DExecutionError, Auto3DUnavailableError, run_auto3d


@dataclass(frozen=True)
class StereoEnergyDecision:
    stereo_id: str
    parent_tautomer_id: str | None
    input_molecule_id: str
    molname: str
    representative_stereo_id: str
    enantiomer_group_id: str
    relationship: str
    energy_kcal_mol: float | None
    relative_energy_kcal_mol: float | None
    energy_rank: int | None
    selected: bool
    reason: str
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class StereoEnergyFilteringResult:
    all_records: list[StereoRecord]
    selected_records: list[StereoRecord]
    rejected_records: list[StereoRecord]
    decisions: list[StereoEnergyDecision]
    collapsed_count: int
    energy_evaluation_count: int


def filter_stereoisomers_with_auto3d(
    records: list[StereoRecord],
    config: RunConfig,
) -> StereoEnergyFilteringResult:
    output_dir = config.output_dir / "stereoisomer_filtering"
    output_dir.mkdir(parents=True, exist_ok=True)
    if not records:
        result = StereoEnergyFilteringResult([], [], [], [], 0, 0)
        write_stereo_energy_outputs(output_dir, result)
        return result

    representative_by_id = _enantiomer_representatives(records, config)
    representatives = _representative_records(records, representative_by_id)
    energies, command, warning = _rank_representatives_with_auto3d(representatives, config, output_dir)
    decisions = _decisions_from_energies(
        records,
        representative_by_id,
        energies,
        command=command,
        fallback_warning=warning,
        config=config,
    )
    selected_ids = {decision.stereo_id for decision in decisions if decision.selected}
    rejected_ids = {decision.stereo_id for decision in decisions if not decision.selected}
    selected = [_annotated_record(record, _decision_for(record.id, decisions)) for record in records if record.id in selected_ids]
    rejected = [_annotated_record(record, _decision_for(record.id, decisions)) for record in records if record.id in rejected_ids]
    annotated_all = selected + rejected
    selected_by_id = {record.id: record for record in selected}
    rejected_by_id = {record.id: record for record in rejected}
    all_records = [selected_by_id.get(record.id) or rejected_by_id.get(record.id) or record for record in records]
    collapsed_count = sum(1 for decision in decisions if decision.relationship == "enantiomer_mapped")
    energy_evaluation_count = len(representatives) if command and not warning else 0
    result = StereoEnergyFilteringResult(
        all_records=all_records,
        selected_records=selected,
        rejected_records=rejected,
        decisions=decisions,
        collapsed_count=collapsed_count,
        energy_evaluation_count=energy_evaluation_count,
    )
    write_stereo_energy_outputs(output_dir, result)
    return result


def write_stereo_energy_outputs(path: Path, result: StereoEnergyFilteringResult) -> None:
    path.mkdir(parents=True, exist_ok=True)
    _write_record_rows(path / "stereoisomers_all.csv", result.all_records, result.decisions)
    _write_record_rows(path / "stereoisomers_selected.csv", result.selected_records, result.decisions)
    _write_record_rows(path / "stereoisomers_rejected.csv", result.rejected_records, result.decisions)
    _write_decision_rows(path / "stereo_energy_ranked.csv", result.decisions)
    _write_group_rows(path / "stereo_enantiomer_groups.csv", result.decisions)
    for name in (
        "stereoisomers_all.csv",
        "stereoisomers_selected.csv",
        "stereoisomers_rejected.csv",
        "stereo_energy_ranked.csv",
        "stereo_enantiomer_groups.csv",
    ):
        _copy_text(path / name, path.parent / name)


def _rank_representatives_with_auto3d(
    records: list[StereoRecord],
    config: RunConfig,
    output_dir: Path,
) -> tuple[dict[str, float], list[str], str | None]:
    if not config.stereoisomer_filtering.enabled:
        return {}, [], "stereoisomer_filtering.enabled=false; retained all stereoisomers"
    if not records:
        return {}, [], None
    try:
        with TemporaryDirectory(prefix="auto3d_stereo_", dir=output_dir) as temp:
            workdir = Path(temp)
            input_smi = workdir / "stereo_representatives.smi"
            input_smi.write_text(
                "".join(f"{record.isomeric_smiles or record.canonical_smiles} {record.id}\n" for record in records),
                encoding="utf-8",
            )
            output_sdf, command = run_auto3d(
                input_smi,
                workdir / "auto3d",
                k=1,
                model=str(config.final_3d.optimizing_engine),
                internal_tautomer_stereo_enum=False,
                max_confs=1,
                patience=config.final_3d.patience,
                use_gpu=config.final_3d.use_gpu,
            )
            energies = _read_auto3d_energies(output_sdf, records)
            if not energies:
                raise Auto3DExecutionError("Auto3D stereo filtering emitted no usable energies")
            return energies, command, None
    except (Auto3DExecutionError, Auto3DUnavailableError) as exc:
        return {}, [], f"Auto3D stereoisomer energy filtering failed; retained all stereoisomers: {exc}"


def _decisions_from_energies(
    records: list[StereoRecord],
    representative_by_id: dict[str, str],
    energies: dict[str, float],
    *,
    command: list[str],
    fallback_warning: str | None,
    config: RunConfig,
) -> list[StereoEnergyDecision]:
    grouped: dict[str | None, list[StereoRecord]] = defaultdict(list)
    for record in records:
        grouped[record.parent_id].append(record)

    selected_ids: set[str] = set()
    rank_by_id: dict[str, int] = {}
    relative_by_id: dict[str, float | None] = {}
    reason_by_id: dict[str, str] = {}

    for group_records in grouped.values():
        representative_ids = sorted({representative_by_id.get(record.id, record.id) for record in group_records})
        ranked_reps = sorted(
            representative_ids,
            key=lambda rid: (float("inf") if energies.get(rid) is None else energies[rid], rid),
        )
        finite = [energies[rid] for rid in ranked_reps if energies.get(rid) is not None]
        minimum = min(finite) if finite else None
        selected_reps: set[str] = set()
        if fallback_warning or minimum is None:
            selected_reps = set(ranked_reps)
        else:
            for rank, representative_id in enumerate(ranked_reps, start=1):
                energy = energies.get(representative_id)
                if energy is None:
                    continue
                relative = energy - minimum
                if (
                    rank <= config.stereoisomer_filtering.keep_top_n_diastereomers
                    and relative <= config.stereoisomer_filtering.stereo_energy_window_kcal_mol
                ):
                    selected_reps.add(representative_id)
        if not selected_reps and ranked_reps:
            selected_reps.add(ranked_reps[0])
        for rank, representative_id in enumerate(ranked_reps, start=1):
            energy = energies.get(representative_id)
            relative = None if energy is None or minimum is None else energy - minimum
            for record in group_records:
                if representative_by_id.get(record.id, record.id) != representative_id:
                    continue
                rank_by_id[record.id] = rank
                relative_by_id[record.id] = relative
                if representative_id in selected_reps:
                    selected_ids.add(record.id)
                    reason_by_id[record.id] = (
                        "selected_without_auto3d_energy" if fallback_warning else "selected_by_auto3d_stereo_energy"
                    )
                elif energy is None:
                    reason_by_id[record.id] = "rejected_missing_auto3d_stereo_energy"
                else:
                    reason_by_id[record.id] = "rejected_by_auto3d_stereo_energy"

    command_warning = f"Auto3D command: {' '.join(str(part) for part in command)}" if command else None
    decisions: list[StereoEnergyDecision] = []
    for record in records:
        representative_id = representative_by_id.get(record.id, record.id)
        relationship = "representative"
        if representative_id != record.id:
            relationship = "enantiomer_mapped"
        warnings = tuple(item for item in (fallback_warning, command_warning) if item)
        decisions.append(
            StereoEnergyDecision(
                stereo_id=record.id,
                parent_tautomer_id=record.parent_id,
                input_molecule_id=record.input_molecule_id,
                molname=record.molname,
                representative_stereo_id=representative_id,
                enantiomer_group_id=_group_id(record, representative_id),
                relationship=relationship,
                energy_kcal_mol=energies.get(representative_id),
                relative_energy_kcal_mol=relative_by_id.get(record.id),
                energy_rank=rank_by_id.get(record.id),
                selected=record.id in selected_ids,
                reason=reason_by_id.get(record.id, "selected_by_default"),
                warnings=warnings,
            )
        )
    return decisions


def _representative_records(
    records: list[StereoRecord],
    representative_by_id: dict[str, str],
) -> list[StereoRecord]:
    by_id = {record.id: record for record in records}
    representative_ids = sorted({representative_by_id.get(record.id, record.id) for record in records})
    return [by_id[record_id] for record_id in representative_ids if record_id in by_id]


def _enantiomer_representatives(records: list[StereoRecord], config: RunConfig) -> dict[str, str]:
    if not (
        config.stereoisomer_filtering.collapse_enantiomers_in_achiral_solvent
        and config.stereoisomer_filtering.run_energy_for_enantiomer_representatives_only
        and not config.stereo_filtering.solvent_is_chiral
    ):
        return {}
    groups: dict[tuple[str | None, int | None, str | None, str], list[StereoRecord]] = defaultdict(list)
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


def _stereo_group_tuple(record: StereoRecord) -> tuple[str | None, int | None, str | None, str]:
    mol = _mol(record)
    achiral_smiles = (
        Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        if mol is not None
        else record.canonical_smiles or ""
    )
    return (record.parent_id, record.formal_charge, record.molecular_formula, achiral_smiles)


def _single_assigned_chiral_center(record: StereoRecord) -> bool:
    centers = _assigned_chiral_centers(record)
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


def _mol(record: StereoRecord) -> Chem.Mol | None:
    if record.rdkit_mol is not None:
        return record.rdkit_mol
    if record.isomeric_smiles:
        return Chem.MolFromSmiles(record.isomeric_smiles)
    if record.canonical_smiles:
        return Chem.MolFromSmiles(record.canonical_smiles)
    return None


def _read_auto3d_energies(output_sdf: Path, records: list[StereoRecord]) -> dict[str, float]:
    by_id = {record.id: record for record in records}
    by_smiles = {record.isomeric_smiles: record for record in records if record.isomeric_smiles}
    energies: dict[str, float] = {}
    supplier = Chem.SDMolSupplier(str(output_sdf), sanitize=True, removeHs=False)
    for molecule in supplier:
        if molecule is None:
            continue
        record = _match_record(molecule, by_id, by_smiles)
        if record is None:
            continue
        energy = _energy_from_mol(molecule)
        if energy is None:
            continue
        previous = energies.get(record.id)
        if previous is None or energy < previous:
            energies[record.id] = energy
    return energies


def _match_record(
    molecule: Chem.Mol,
    by_id: dict[str, StereoRecord],
    by_smiles: dict[str, StereoRecord],
) -> StereoRecord | None:
    for key in ("DSVR_STEREO_ID", "stereo_id", "ID", "_Name"):
        if molecule.HasProp(key):
            value = molecule.GetProp(key).strip()
            if value in by_id:
                return by_id[value]
            token = value.split()[0] if value else ""
            if token in by_id:
                return by_id[token]
    smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), canonical=True, isomericSmiles=True)
    return by_smiles.get(smiles)


def _energy_from_mol(molecule: Chem.Mol) -> float | None:
    for key in (
        "E_kcal_mol",
        "energy_kcal_mol",
        "E_stereo_relative(kcal/mol)",
        "E",
        "Energy",
        "energy",
        "Auto3D_energy",
    ):
        if molecule.HasProp(key):
            try:
                return float(molecule.GetProp(key))
            except ValueError:
                continue
    return None


def _annotated_record(record: StereoRecord, decision: StereoEnergyDecision) -> StereoRecord:
    metadata = dict(record.metadata)
    metadata["stereo_energy_filtering"] = {
        "representative_stereo_id": decision.representative_stereo_id,
        "enantiomer_group_id": decision.enantiomer_group_id,
        "relationship": decision.relationship,
        "energy_kcal_mol": decision.energy_kcal_mol,
        "relative_energy_kcal_mol": decision.relative_energy_kcal_mol,
        "energy_rank": decision.energy_rank,
        "selected": decision.selected,
        "reason": decision.reason,
    }
    warnings = sorted({*record.warnings, *decision.warnings})
    if decision.relationship == "enantiomer_mapped":
        warnings.append(
            "Auto3D stereo energy was mapped from an enantiomeric representative in achiral solvent."
        )
    return record.model_copy(update={"metadata": metadata, "warnings": warnings})


def _decision_for(stereo_id: str, decisions: list[StereoEnergyDecision]) -> StereoEnergyDecision:
    for decision in decisions:
        if decision.stereo_id == stereo_id:
            return decision
    raise KeyError(stereo_id)


def _group_id(record: StereoRecord, representative_id: str) -> str:
    return f"{record.parent_id or record.input_molecule_id}:{representative_id}"


def _write_record_rows(
    path: Path,
    records: list[StereoRecord],
    decisions: list[StereoEnergyDecision],
) -> None:
    by_id = {decision.stereo_id: decision for decision in decisions}
    columns = [
        "id",
        "parent_id",
        "input_molecule_id",
        "molname",
        "canonical_smiles",
        "isomeric_smiles",
        "molecular_formula",
        "formal_charge",
        "explicit_proton_count",
        "representative_stereo_id",
        "enantiomer_group_id",
        "relationship",
        "energy_kcal_mol",
        "relative_energy_kcal_mol",
        "energy_rank",
        "selected",
        "reason",
        "warnings",
    ]
    rows: list[dict[str, Any]] = []
    for record in records:
        decision = by_id.get(record.id)
        rows.append(
            {
                "id": record.id,
                "parent_id": record.parent_id,
                "input_molecule_id": record.input_molecule_id,
                "molname": record.molname,
                "canonical_smiles": record.canonical_smiles,
                "isomeric_smiles": record.isomeric_smiles,
                "molecular_formula": record.molecular_formula,
                "formal_charge": record.formal_charge,
                "explicit_proton_count": record.explicit_proton_count,
                "representative_stereo_id": decision.representative_stereo_id if decision else record.id,
                "enantiomer_group_id": decision.enantiomer_group_id if decision else "",
                "relationship": decision.relationship if decision else "",
                "energy_kcal_mol": decision.energy_kcal_mol if decision else None,
                "relative_energy_kcal_mol": decision.relative_energy_kcal_mol if decision else None,
                "energy_rank": decision.energy_rank if decision else None,
                "selected": decision.selected if decision else True,
                "reason": decision.reason if decision else "",
                "warnings": " | ".join(record.warnings),
            }
        )
    _write_rows(path, columns, rows)


def _write_decision_rows(path: Path, decisions: list[StereoEnergyDecision]) -> None:
    columns = [
        "stereo_id",
        "parent_tautomer_id",
        "input_molecule_id",
        "molname",
        "representative_stereo_id",
        "enantiomer_group_id",
        "relationship",
        "energy_kcal_mol",
        "relative_energy_kcal_mol",
        "energy_rank",
        "selected",
        "reason",
        "warnings",
    ]
    _write_rows(
        path,
        columns,
        [
            {
                **decision.__dict__,
                "warnings": " | ".join(decision.warnings),
            }
            for decision in decisions
        ],
    )


def _write_group_rows(path: Path, decisions: list[StereoEnergyDecision]) -> None:
    columns = [
        "enantiomer_group_id",
        "representative_stereo_id",
        "member_stereo_ids",
        "energy_evaluated_once",
        "member_count",
    ]
    groups: dict[str, list[StereoEnergyDecision]] = defaultdict(list)
    for decision in decisions:
        groups[decision.enantiomer_group_id].append(decision)
    rows = []
    for group_id, members in sorted(groups.items()):
        representative_id = sorted({member.representative_stereo_id for member in members})[0]
        rows.append(
            {
                "enantiomer_group_id": group_id,
                "representative_stereo_id": representative_id,
                "member_stereo_ids": " | ".join(member.stereo_id for member in members),
                "energy_evaluated_once": any(member.relationship == "enantiomer_mapped" for member in members),
                "member_count": len(members),
            }
        )
    _write_rows(path, columns, rows)


def _write_rows(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column) for column in columns})


def _copy_text(source: Path, target: Path) -> None:
    if source.resolve() == target.resolve():
        return
    target.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
