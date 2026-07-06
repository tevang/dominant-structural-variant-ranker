from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.config import RunConfig
from dsvr.models import MoleculeInput, ProtomerRecord, make_protomer_id
from dsvr.runners.molscrub_runner import generate_molscrub_candidates


@dataclass(frozen=True)
class ProtomerPlausibility:
    score: float
    reasons: list[str]
    warnings: list[str]
    selected: bool
    selection_reason: str | None


def describe_protonation_scope(ph: float) -> str:
    return f"candidate-generation pH={ph:g}; no rigorous pH population correction"


def generate_protomer_candidates(
    mol_record: MoleculeInput,
    config: RunConfig,
) -> list[ProtomerRecord]:
    if not config.protonation.enabled:
        output_dir = config.output_dir / "enumeration" / "protomers"
        output_dir.mkdir(parents=True, exist_ok=True)
        records = _records_from_candidates(
            mol_record,
            [Chem.Mol(mol_record.rdkit_mol)],
            config=config,
            source_software="input-fallback",
            source_command="protonation.enabled=false",
            output_dir=output_dir,
            fallback_warning="protonation.enabled=false; retained input state only",
        )
        return records

    ph_low = config.chemistry.ph_low if config.chemistry.ph_low is not None else config.chemistry.ph
    ph_high = config.chemistry.ph_high if config.chemistry.ph_high is not None else config.chemistry.ph
    raw_candidates, source_software, source_command = generate_molscrub_candidates(
        mol_record.rdkit_mol,
        ph_low=ph_low,
        ph_high=ph_high,
        skip_gen3d=config.protonation.skip_gen3d_in_molscrub,
        timeout_seconds=config.protonation.timeout_seconds_per_molecule,
    )
    output_dir = config.output_dir / "enumeration" / "protomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    return _records_from_candidates(
        mol_record,
        raw_candidates,
        config=config,
        source_software=source_software,
        source_command=source_command,
        output_dir=output_dir,
    )


def _records_from_candidates(
    mol_record: MoleculeInput,
    candidates: list[Chem.Mol],
    *,
    config: RunConfig,
    source_software: str,
    source_command: str,
    output_dir: Path,
    fallback_warning: str | None = None,
) -> list[ProtomerRecord]:
    if not candidates:
        fallback_warning = fallback_warning or "molscrub returned no valid state; retained input molecule"
    unique_candidates, duplicate_rows = _dedupe_candidates(mol_record, candidates, config)
    if not unique_candidates:
        unique_candidates = [Chem.Mol(mol_record.rdkit_mol)]
        fallback_warning = fallback_warning or "molscrub returned no valid state; retained input molecule"

    selected, rejected_rows, plausibility_by_key = _select_plausible_protomers(
        mol_record,
        unique_candidates,
        config,
        duplicate_rows,
    )
    if not selected:
        fallback = Chem.Mol(mol_record.rdkit_mol)
        key = _candidate_key(fallback)
        selected = [fallback]
        plausibility_by_key[key] = ProtomerPlausibility(
            score=0.0,
            reasons=["input fallback retained"],
            warnings=[fallback_warning or "no protomer selected; retained input molecule"],
            selected=True,
            selection_reason="fallback_input_state",
        )

    records: list[ProtomerRecord] = []
    for index, candidate in enumerate(selected, start=1):
        canonical_smiles = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True)
        formula = _formula(candidate)
        charge = Chem.GetFormalCharge(candidate)
        proton_count = _explicit_proton_count(candidate)
        plausibility = plausibility_by_key[_candidate_key(candidate)]
        warnings = [
            "molscrub pH influence is candidate generation/filtering only; no rigorous pH "
            "population prediction is implied."
        ]
        warnings.extend(plausibility.warnings)
        if fallback_warning and fallback_warning not in warnings:
            warnings.append(fallback_warning)
        metadata = {
            "ph_low": config.chemistry.ph_low,
            "ph_high": config.chemistry.ph_high,
            "target_ph": config.chemistry.ph,
            "solvent": config.chemistry.solvent,
            "candidate_generation_only": True,
            "score_is_population_estimate": False,
            "dedupe_key": {
                "formula": formula,
                "formal_charge": charge,
                "canonical_smiles": canonical_smiles,
                "isomeric_smiles": isomeric_smiles,
            },
            "plausibility": asdict(plausibility),
        }
        records.append(
            ProtomerRecord(
                id=make_protomer_id(
                    mol_record.input_id,
                    index,
                    canonical_smiles,
                    isomeric_smiles,
                    metadata,
                ),
                parent_id=mol_record.input_id,
                input_molecule_id=mol_record.input_id,
                molname=mol_record.molname,
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                molecular_formula=formula,
                formal_charge=charge,
                explicit_proton_count=proton_count,
                source_software=source_software,
                source_command=source_command,
                source_python_function="dsvr.chemistry.protonation.generate_protomer_candidates",
                output_paths=[
                    output_dir / f"{mol_record.input_id}_protomers.sdf",
                    output_dir / f"{mol_record.input_id}_protomers.csv",
                    output_dir / "protomers_selected.csv",
                    output_dir / "protomers_rejected.csv",
                    output_dir / "protonation_warnings.jsonl",
                ],
                warnings=warnings,
                metadata=metadata,
                protomer_index=index,
                rdkit_mol=candidate,
            )
        )

    _write_protomer_sdf(output_dir / f"{mol_record.input_id}_protomers.sdf", records)
    _write_protomer_csv(output_dir / f"{mol_record.input_id}_protomers.csv", records)
    _write_audit_outputs(output_dir, mol_record, unique_candidates, records, rejected_rows, plausibility_by_key)
    return records


def _dedupe_candidates(
    mol_record: MoleculeInput,
    candidates: list[Chem.Mol],
    config: RunConfig,
) -> tuple[list[Chem.Mol], list[dict[str, Any]]]:
    candidate_list = list(candidates)
    if config.protonation.keep_input_state and mol_record.rdkit_mol is not None:
        candidate_list.insert(0, Chem.Mol(mol_record.rdkit_mol))

    seen: set[tuple[str, int, str, str]] = set()
    unique: list[Chem.Mol] = []
    duplicate_rows: list[dict[str, Any]] = []
    for index, candidate in enumerate(candidate_list, start=1):
        if candidate is None:
            duplicate_rows.append({"input_index": index, "selected": False, "reason": "invalid_molecule"})
            continue
        try:
            key = _candidate_key(candidate)
        except Exception as exc:
            duplicate_rows.append(
                {"input_index": index, "selected": False, "reason": f"invalid_molecule: {exc}"}
            )
            continue
        if key in seen:
            duplicate_rows.append(
                {
                    "input_index": index,
                    "selected": False,
                    "reason": "duplicate_dedupe_key",
                    "dedupe_key": _json_key(key),
                }
            )
            continue
        seen.add(key)
        unique.append(candidate)
    return unique, duplicate_rows


def _select_plausible_protomers(
    mol_record: MoleculeInput,
    candidates: list[Chem.Mol],
    config: RunConfig,
    duplicate_rows: list[dict[str, Any]],
) -> tuple[list[Chem.Mol], list[dict[str, Any]], dict[tuple[str, int, str, str], ProtomerPlausibility]]:
    cap = config.protonation.max_protomers_per_molecule
    scored = [(candidate, _score_candidate(candidate, mol_record, config)) for candidate in candidates]
    selected_keys: set[tuple[str, int, str, str]] = set()
    selected: list[Chem.Mol] = []

    if config.protonation.keep_input_state:
        input_key = _candidate_key(mol_record.rdkit_mol)
        for candidate, plausibility in scored:
            if _candidate_key(candidate) == input_key and len(selected) < cap:
                selected.append(candidate)
                selected_keys.add(input_key)
                scored = [
                    (mol, _select_plausibility(pl, "input_reference_state"))
                    if _candidate_key(mol) == input_key
                    else (mol, pl)
                    for mol, pl in scored
                ]
                break

    if config.protonation.keep_best_per_charge:
        best_by_charge: dict[int, tuple[Chem.Mol, ProtomerPlausibility]] = {}
        for candidate, plausibility in scored:
            charge = Chem.GetFormalCharge(candidate)
            current = best_by_charge.get(charge)
            if current is None or plausibility.score < current[1].score:
                best_by_charge[charge] = (candidate, plausibility)
        for candidate, plausibility in sorted(
            best_by_charge.values(), key=lambda item: (item[1].score, abs(Chem.GetFormalCharge(item[0])))
        ):
            key = _candidate_key(candidate)
            if key not in selected_keys and len(selected) < cap:
                selected.append(candidate)
                selected_keys.add(key)
                scored = [
                    (mol, _select_plausibility(pl, "best_representative_per_charge"))
                    if _candidate_key(mol) == key
                    else (mol, pl)
                    for mol, pl in scored
                ]

    for candidate, plausibility in sorted(scored, key=lambda item: (item[1].score, item[0].GetNumAtoms())):
        key = _candidate_key(candidate)
        if key not in selected_keys and len(selected) < cap:
            selected.append(candidate)
            selected_keys.add(key)
            scored = [
                (mol, _select_plausibility(pl, "score_ranked_fill")) if _candidate_key(mol) == key else (mol, pl)
                for mol, pl in scored
            ]

    plausibility_by_key: dict[tuple[str, int, str, str], ProtomerPlausibility] = {}
    rejected_rows = list(duplicate_rows)
    for index, (candidate, plausibility) in enumerate(scored, start=1):
        key = _candidate_key(candidate)
        if key in selected_keys:
            plausibility_by_key[key] = plausibility if plausibility.selected else _select_plausibility(plausibility, "selected")
            continue
        reason = "beyond_max_protomers_per_molecule"
        same_charge_selected = any(Chem.GetFormalCharge(item) == Chem.GetFormalCharge(candidate) for item in selected)
        if same_charge_selected:
            reason = "lower_scoring_same_charge_state"
        rejected = ProtomerPlausibility(
            score=plausibility.score,
            reasons=plausibility.reasons,
            warnings=plausibility.warnings,
            selected=False,
            selection_reason=reason,
        )
        plausibility_by_key[key] = rejected
        rejected_rows.append(_candidate_row(candidate, mol_record, index, rejected, reason))
    return selected, rejected_rows, plausibility_by_key


def _score_candidate(candidate: Chem.Mol, mol_record: MoleculeInput, config: RunConfig) -> ProtomerPlausibility:
    score = 0.0
    reasons: list[str] = []
    warnings: list[str] = []
    charge = Chem.GetFormalCharge(candidate)
    if abs(charge) > 0:
        penalty = 1.5 * abs(charge)
        score += penalty
        reasons.append(f"formal_charge_penalty={penalty:g}")
    if abs(charge) > 1:
        penalty = 2.0 * (abs(charge) - 1)
        score += penalty
        reasons.append(f"excessive_formal_charge_penalty={penalty:g}")
    if _is_nonpolar_solvent(config.chemistry.solvent) and charge != 0:
        score += 2.0
        reasons.append("charged_state_penalty_in_nonpolar_solvent=2")
    if _candidate_key(candidate) == _candidate_key(mol_record.rdkit_mol):
        score -= 1.0
        reasons.append("input_reference_state_bonus=-1")
    score_prop = _first_float_prop(candidate, ["molscrub_score", "MOLSCRUB_SCORE", "score", "probability"])
    if score_prop is not None:
        score -= max(0.0, min(score_prop, 1.0))
        reasons.append("molscrub_score_used_as_priority_hint")
    else:
        warnings.append("molscrub did not expose reliable pH abundance/probability; no pKa score invented")
    return ProtomerPlausibility(
        score=score,
        reasons=reasons or ["neutral_score"],
        warnings=warnings,
        selected=False,
        selection_reason=None,
    )


def _select_plausibility(plausibility: ProtomerPlausibility, reason: str) -> ProtomerPlausibility:
    return ProtomerPlausibility(
        score=plausibility.score,
        reasons=plausibility.reasons,
        warnings=plausibility.warnings,
        selected=True,
        selection_reason=reason,
    )


def _write_audit_outputs(
    output_dir: Path,
    mol_record: MoleculeInput,
    all_candidates: list[Chem.Mol],
    selected_records: list[ProtomerRecord],
    rejected_rows: list[dict[str, Any]],
    plausibility_by_key: dict[tuple[str, int, str, str], ProtomerPlausibility],
) -> None:
    all_rows = []
    for index, candidate in enumerate(all_candidates, start=1):
        key = _candidate_key(candidate)
        plausibility = plausibility_by_key.get(key) or _score_candidate(candidate, mol_record, RunConfig())
        all_rows.append(_candidate_row(candidate, mol_record, index, plausibility, plausibility.selection_reason))
    _append_csv(output_dir / "protomers_all.csv", all_rows)
    _append_csv(
        output_dir / "protomers_selected.csv",
        [_record_row(record) for record in selected_records],
    )
    _append_csv(output_dir / "protomers_rejected.csv", rejected_rows)
    warnings_path = output_dir / "protonation_warnings.jsonl"
    with warnings_path.open("a", encoding="utf-8") as handle:
        for record in selected_records:
            for warning in record.warnings:
                handle.write(
                    json.dumps(
                        {
                            "input_molecule_id": record.input_molecule_id,
                            "protomer_id": record.id,
                            "warning": warning,
                        },
                        sort_keys=True,
                    )
                    + "\n"
                )


def _append_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        if not path.exists():
            path.write_text("", encoding="utf-8")
        return
    columns = sorted({key for row in rows for key in row})
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def _candidate_row(
    candidate: Chem.Mol,
    mol_record: MoleculeInput,
    index: int,
    plausibility: ProtomerPlausibility,
    reason: str | None,
) -> dict[str, Any]:
    canonical_smiles = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=False)
    isomeric_smiles = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True)
    return {
        "input_molecule_id": mol_record.input_id,
        "molname": mol_record.molname,
        "candidate_index": index,
        "canonical_smiles": canonical_smiles,
        "isomeric_smiles": isomeric_smiles,
        "molecular_formula": _formula(candidate),
        "formal_charge": Chem.GetFormalCharge(candidate),
        "score": plausibility.score,
        "selected": plausibility.selected,
        "selection_reason": plausibility.selection_reason or reason or "",
        "reasons": " | ".join(plausibility.reasons),
        "warnings": " | ".join(plausibility.warnings),
        "score_is_population_estimate": False,
    }


def _record_row(record: ProtomerRecord) -> dict[str, Any]:
    plausibility = record.metadata.get("plausibility", {})
    return {
        "input_molecule_id": record.input_molecule_id,
        "protomer_id": record.id,
        "molname": record.molname,
        "canonical_smiles": record.canonical_smiles,
        "isomeric_smiles": record.isomeric_smiles,
        "molecular_formula": record.molecular_formula,
        "formal_charge": record.formal_charge,
        "score": plausibility.get("score"),
        "selected": True,
        "selection_reason": plausibility.get("selection_reason"),
        "warnings": " | ".join(record.warnings),
        "score_is_population_estimate": False,
    }


def _write_protomer_sdf(path: Path, records: list[ProtomerRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        mol = Chem.Mol(record.rdkit_mol)
        mol.SetProp("_Name", record.id)
        for key, value in {
            "DSVR_STAGE": record.stage_name,
            "DSVR_INPUT_ID": record.input_molecule_id,
            "DSVR_PARENT_ID": record.parent_id or "",
            "DSVR_PROTOMER_ID": record.id,
            "DSVR_MOLNAME": record.molname,
            "DSVR_CANONICAL_SMILES": record.canonical_smiles or "",
            "DSVR_ISOMERIC_SMILES": record.isomeric_smiles or "",
            "DSVR_FORMULA": record.molecular_formula or "",
            "DSVR_FORMAL_CHARGE": str(record.formal_charge),
            "DSVR_EXPLICIT_PROTON_COUNT": str(record.explicit_proton_count),
            "DSVR_PH_SCOPE": "candidate_generation_filtering_only_not_population",
            "DSVR_PROTOMER_SCORE": str(record.metadata.get("plausibility", {}).get("score", "")),
        }.items():
            mol.SetProp(key, value)
        writer.write(mol)
    writer.close()


def _write_protomer_csv(path: Path, records: list[ProtomerRecord]) -> None:
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
        "stage_name",
        "source_software",
        "source_command",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = record.model_dump(mode="json")
            row["warnings"] = " | ".join(record.warnings)
            writer.writerow({column: row.get(column) for column in columns})


def _candidate_key(molecule: Chem.Mol) -> tuple[str, int, str, str]:
    return (
        _formula(molecule),
        Chem.GetFormalCharge(molecule),
        Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False),
        Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
    )


def _json_key(key: tuple[str, int, str, str]) -> str:
    return json.dumps(
        {"formula": key[0], "formal_charge": key[1], "canonical_smiles": key[2], "isomeric_smiles": key[3]},
        sort_keys=True,
    )


def _formula(molecule: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(molecule)


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)


def _is_nonpolar_solvent(solvent: str) -> bool:
    return solvent.strip().lower() in {"benzene", "chloroform", "ether", "hexane", "toluene"}


def _first_float_prop(molecule: Chem.Mol, names: list[str]) -> float | None:
    for name in names:
        if not molecule.HasProp(name):
            continue
        try:
            return float(molecule.GetProp(name))
        except ValueError:
            continue
    return None
