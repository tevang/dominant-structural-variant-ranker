from __future__ import annotations

import csv
import multiprocessing as mp
import queue
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from dsvr.config import RunConfig
from dsvr.models import ProtomerRecord, TautomerRecord, make_tautomer_id
from dsvr.runners.auto3d_runner import Auto3DExecutionError, Auto3DUnavailableError, run_auto3d


class Auto3DTautomerFilteringError(RuntimeError):
    """Raised when Auto3D tautomer triage cannot produce ranked candidates."""


class RdkitTautomerFilteringTimeout(RuntimeError):
    """Raised when RDKit tautomer candidate enumeration times out."""


@dataclass(frozen=True)
class _Candidate:
    index: int
    tautomer_id: str
    molecule: Chem.Mol
    canonical_smiles: str
    isomeric_smiles: str
    rdkit_score: float | None
    is_input_tautomer: bool
    is_canonical_tautomer: bool


@dataclass(frozen=True)
class _RankedCandidate:
    candidate: _Candidate
    energy_kcal_mol: float | None
    relative_energy_kcal_mol: float | None
    auto3d_rank: int | None
    selected: bool
    reason: str
    source: str
    warnings: tuple[str, ...] = ()


def filter_tautomers_with_auto3d(
    protomer_records: list[ProtomerRecord],
    config: RunConfig,
) -> list[TautomerRecord]:
    output_dir = config.output_dir / "enumeration" / "tautomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    selected_records: list[TautomerRecord] = []
    all_rows: list[dict[str, Any]] = []
    ranked_rows: list[dict[str, Any]] = []
    selected_rows: list[dict[str, Any]] = []
    rejected_rows: list[dict[str, Any]] = []

    for protomer in protomer_records:
        candidates, candidate_warning = _enumerate_candidates(protomer, config)
        ranked = _rank_or_fallback(protomer, candidates, config, output_dir, candidate_warning)
        records = _records_from_ranked(
            protomer,
            [item for item in ranked if item.selected],
            config,
            output_dir,
        )
        selected_records.extend(records)
        _write_tautomer_sdf(output_dir / f"{protomer.id}_tautomers.sdf", records)
        _write_tautomer_csv(output_dir / f"{protomer.id}_tautomers.csv", records)
        all_rows.extend(_candidate_row(protomer, item, candidate_warning) for item in candidates)
        ranked_rows.extend(_ranked_row(protomer, item) for item in ranked)
        selected_rows.extend(_ranked_row(protomer, item) for item in ranked if item.selected)
        rejected_rows.extend(_ranked_row(protomer, item) for item in ranked if not item.selected)

    _write_rows(output_dir / "tautomers_all_pre_auto3d.csv", all_rows)
    _write_rows(output_dir / "tautomers_auto3d_ranked.csv", ranked_rows)
    _write_rows(output_dir / "tautomers_selected.csv", selected_rows)
    _write_rows(output_dir / "tautomers_rejected.csv", rejected_rows)
    return selected_records


def _enumerate_candidates(
    protomer: ProtomerRecord,
    config: RunConfig,
) -> tuple[list[_Candidate], str | None]:
    warning: str | None = None
    input_mol = Chem.Mol(protomer.rdkit_mol)
    params = rdMolStandardize.CleanupParameters()
    params.maxTautomers = config.tautomer_filtering.max_rdkit_tautomers_before_auto3d
    params.maxTransforms = config.enumeration.max_tautomer_transforms
    params.tautomerRemoveBondStereo = config.enumeration.tautomer_remove_bond_stereo
    params.tautomerRemoveSp3Stereo = config.enumeration.tautomer_remove_sp3_stereo
    params.tautomerReassignStereo = config.enumeration.tautomer_reassign_stereo
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    try:
        molblocks = _enumerate_molblocks_with_timeout(
            input_mol,
            timeout_seconds=config.tautomer_filtering.rdkit_tautomer_timeout_seconds,
            max_tautomers=config.tautomer_filtering.max_rdkit_tautomers_before_auto3d,
            max_transforms=config.enumeration.max_tautomer_transforms,
            remove_bond_stereo=config.enumeration.tautomer_remove_bond_stereo,
            remove_sp3_stereo=config.enumeration.tautomer_remove_sp3_stereo,
            reassign_stereo=config.enumeration.tautomer_reassign_stereo,
        )
        molecules = [
            molecule
            for molecule in (Chem.MolFromMolBlock(block, sanitize=True, removeHs=False) for block in molblocks)
            if molecule is not None
        ]
    except RdkitTautomerFilteringTimeout:
        molecules = [input_mol]
        warning = "RDKit tautomer enumeration timeout; retained input tautomer only"
    except RuntimeError as exc:
        molecules = [input_mol]
        warning = f"RDKit tautomer enumeration failed; retained input tautomer only: {exc}"

    if config.tautomer_filtering.keep_input_tautomer:
        molecules.append(input_mol)
    canonical_mol = enumerator.Canonicalize(input_mol)
    if canonical_mol is not None:
        molecules.append(canonical_mol)
    return _dedupe_candidates(protomer, molecules, enumerator), warning


def _enumerate_molblocks_with_timeout(
    molecule: Chem.Mol,
    *,
    timeout_seconds: int,
    max_tautomers: int,
    max_transforms: int,
    remove_bond_stereo: bool,
    remove_sp3_stereo: bool,
    reassign_stereo: bool,
) -> list[str]:
    output_queue: mp.Queue = mp.Queue(maxsize=1)
    process = mp.Process(
        target=_tautomer_worker,
        args=(
            Chem.MolToMolBlock(molecule),
            max_tautomers,
            max_transforms,
            remove_bond_stereo,
            remove_sp3_stereo,
            reassign_stereo,
            output_queue,
        ),
    )
    process.start()
    process.join(timeout_seconds)
    if process.is_alive():
        process.terminate()
        process.join(2)
        if process.is_alive():
            process.kill()
            process.join(2)
        raise RdkitTautomerFilteringTimeout("RDKit tautomer enumeration timed out")
    try:
        payload = output_queue.get_nowait()
    except queue.Empty as exc:
        raise RuntimeError("RDKit tautomer worker produced no output") from exc
    if payload.get("status") != "ok":
        raise RuntimeError(str(payload.get("error", "unknown tautomer worker error")))
    return list(payload.get("molblocks", []))


def _tautomer_worker(
    molblock: str,
    max_tautomers: int,
    max_transforms: int,
    remove_bond_stereo: bool,
    remove_sp3_stereo: bool,
    reassign_stereo: bool,
    output_queue: mp.Queue,
) -> None:
    try:
        molecule = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
        if molecule is None:
            raise ValueError("could not parse tautomer worker molecule")
        params = rdMolStandardize.CleanupParameters()
        params.maxTautomers = max_tautomers
        params.maxTransforms = max_transforms
        params.tautomerRemoveBondStereo = remove_bond_stereo
        params.tautomerRemoveSp3Stereo = remove_sp3_stereo
        params.tautomerReassignStereo = reassign_stereo
        enumerator = rdMolStandardize.TautomerEnumerator(params)
        tautomers = list(enumerator.Enumerate(molecule))[:max_tautomers]
        output_queue.put({"status": "ok", "molblocks": [Chem.MolToMolBlock(item) for item in tautomers]})
    except Exception as exc:  # pragma: no cover - exercised through parent process.
        output_queue.put({"status": "error", "error": str(exc)})


def _dedupe_candidates(
    protomer: ProtomerRecord,
    molecules: list[Chem.Mol],
    enumerator: rdMolStandardize.TautomerEnumerator,
) -> list[_Candidate]:
    input_isomeric = Chem.MolToSmiles(protomer.rdkit_mol, canonical=True, isomericSmiles=True)
    canonical_mol = enumerator.Canonicalize(Chem.Mol(protomer.rdkit_mol))
    canonical_isomeric = (
        Chem.MolToSmiles(canonical_mol, canonical=True, isomericSmiles=True)
        if canonical_mol is not None
        else input_isomeric
    )
    seen: set[tuple[str, str]] = set()
    candidates: list[_Candidate] = []
    for molecule in molecules:
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        key = (canonical_smiles, isomeric_smiles)
        if key in seen:
            continue
        seen.add(key)
        index = len(candidates) + 1
        metadata = {"auto3d_tautomer_filtering": True, "tautomer_filtering_stage": "pre_stereo"}
        candidates.append(
            _Candidate(
                index=index,
                tautomer_id=make_tautomer_id(protomer.id, index, canonical_smiles, isomeric_smiles, metadata),
                molecule=molecule,
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                rdkit_score=_score_tautomer(enumerator, molecule),
                is_input_tautomer=isomeric_smiles == input_isomeric,
                is_canonical_tautomer=isomeric_smiles == canonical_isomeric,
            )
        )
    return candidates


def _rank_or_fallback(
    protomer: ProtomerRecord,
    candidates: list[_Candidate],
    config: RunConfig,
    output_dir: Path,
    candidate_warning: str | None,
) -> list[_RankedCandidate]:
    if not candidates:
        return []
    if len(candidates) == 1:
        return _single_candidate_rank(candidates[0], candidate_warning or "single RDKit tautomer candidate")
    if not config.tautomer_filtering.enabled:
        return _fallback_rank(candidates, "tautomer_filtering.enabled=false", config)
    try:
        return _rank_with_auto3d(protomer, candidates, config, output_dir)
    except (Auto3DExecutionError, Auto3DUnavailableError, Auto3DTautomerFilteringError) as exc:
        reason = f"Auto3D tautomer filtering failed; RDKit fallback used: {exc}"
        if candidate_warning:
            reason = f"{candidate_warning}; {reason}"
        return _fallback_rank(candidates, reason, config)


def _rank_with_auto3d(
    protomer: ProtomerRecord,
    candidates: list[_Candidate],
    config: RunConfig,
    output_dir: Path,
) -> list[_RankedCandidate]:
    with TemporaryDirectory(prefix=f"{protomer.id}_auto3d_tautomers_", dir=output_dir) as temp:
        workdir = Path(temp)
        input_smi = workdir / "tautomer_candidates.smi"
        input_smi.write_text(
            "".join(f"{candidate.isomeric_smiles} {candidate.tautomer_id}\n" for candidate in candidates),
            encoding="utf-8",
        )
        output_sdf, command = _run_auto3d_for_tautomers(input_smi, workdir, config)
        energies = _read_auto3d_energies(output_sdf, candidates)
        if not energies:
            raise Auto3DTautomerFilteringError("Auto3D output SDF did not contain usable tautomer energies")
    ranked = sorted(
        ((energy, candidate.isomeric_smiles, candidate) for candidate in candidates if (energy := energies.get(candidate.tautomer_id)) is not None),
        key=lambda item: (item[0], item[1]),
    )
    if not ranked:
        raise Auto3DTautomerFilteringError("Auto3D did not emit energies for any RDKit tautomer candidate")
    minimum = ranked[0][0]
    selected_ids = _selected_candidate_ids(ranked, config)
    selected_ids = _rescue_input_tautomer(selected_ids, ranked, candidates, config)
    rank_by_id = {candidate.tautomer_id: rank for rank, (_energy, _smiles, candidate) in enumerate(ranked, start=1)}
    command_text = " ".join(str(part) for part in command)
    results: list[_RankedCandidate] = []
    for candidate in candidates:
        energy = energies.get(candidate.tautomer_id)
        relative = None if energy is None else energy - minimum
        selected = candidate.tautomer_id in selected_ids
        if energy is None:
            reason = "rejected_missing_auto3d_energy"
        elif selected:
            reason = "selected_by_auto3d_energy_filter"
        else:
            reason = "rejected_by_auto3d_energy_filter"
        results.append(
            _RankedCandidate(
                candidate=candidate,
                energy_kcal_mol=energy,
                relative_energy_kcal_mol=relative,
                auto3d_rank=rank_by_id.get(candidate.tautomer_id),
                selected=selected,
                reason=reason,
                source="auto3d",
                warnings=(f"Auto3D command: {command_text}",),
            )
        )
    return results


def _run_auto3d_for_tautomers(input_smi: Path, workdir: Path, config: RunConfig) -> tuple[Path, list[str]]:
    errors: list[str] = []
    engines = [
        config.tautomer_filtering.optimizing_engine,
        config.tautomer_filtering.fallback_optimizing_engine,
        "AIMNet2",
        "AIMNET",
    ]
    seen: set[str] = set()
    for engine in engines:
        engine_name = str(engine)
        if engine_name in seen:
            continue
        seen.add(engine_name)
        try:
            return run_auto3d(
                input_smi,
                workdir / f"auto3d_{engine_name}",
                k=config.tautomer_filtering.auto3d_max_confs_per_tautomer,
                model=engine_name,
                internal_tautomer_stereo_enum=False,
                max_confs=config.tautomer_filtering.auto3d_max_confs_per_tautomer,
                patience=config.tautomer_filtering.auto3d_patience,
                use_gpu=config.tautomer_filtering.use_gpu,
            )
        except (Auto3DExecutionError, Auto3DUnavailableError) as exc:
            errors.append(f"{engine_name}: {exc}")
    raise Auto3DExecutionError("Auto3D tautomer ranking failed. Tried engines:\n" + "\n".join(errors))


def _read_auto3d_energies(output_sdf: Path, candidates: list[_Candidate]) -> dict[str, float]:
    by_id = {candidate.tautomer_id: candidate for candidate in candidates}
    by_smiles = {candidate.isomeric_smiles: candidate for candidate in candidates}
    energies: dict[str, float] = {}
    supplier = Chem.SDMolSupplier(str(output_sdf), sanitize=True, removeHs=False)
    for molecule in supplier:
        if molecule is None:
            continue
        candidate = _match_candidate(molecule, by_id, by_smiles)
        if candidate is None:
            continue
        energy = _energy_from_mol(molecule)
        if energy is None:
            continue
        previous = energies.get(candidate.tautomer_id)
        if previous is None or energy < previous:
            energies[candidate.tautomer_id] = energy
    return energies


def _match_candidate(molecule: Chem.Mol, by_id: dict[str, _Candidate], by_smiles: dict[str, _Candidate]) -> _Candidate | None:
    for key in ("DSVR_TAUTOMER_ID", "tautomer_id", "ID", "_Name"):
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
        "E_tautomer_relative(kcal/mol)",
        "E_kcal_mol",
        "energy_kcal_mol",
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


def _selected_candidate_ids(ranked: list[tuple[float, str, _Candidate]], config: RunConfig) -> set[str]:
    selected: set[str] = set()
    minimum = ranked[0][0]
    for rank, (energy, _smiles, candidate) in enumerate(ranked, start=1):
        if rank > config.tautomer_filtering.tauto_k:
            continue
        if energy - minimum <= config.tautomer_filtering.tauto_window_kcal_mol:
            selected.add(candidate.tautomer_id)
    if not selected:
        selected.add(ranked[0][2].tautomer_id)
    return selected


def _rescue_input_tautomer(
    selected_ids: set[str],
    ranked: list[tuple[float, str, _Candidate]],
    candidates: list[_Candidate],
    config: RunConfig,
) -> set[str]:
    if not config.tautomer_filtering.keep_input_tautomer:
        return selected_ids
    input_candidate = next((candidate for candidate in candidates if candidate.is_input_tautomer), None)
    if input_candidate is None or input_candidate.tautomer_id in selected_ids:
        return selected_ids
    selected = set(selected_ids)
    if len(selected) >= config.tautomer_filtering.tauto_k:
        selected_ranked = [item for item in ranked if item[2].tautomer_id in selected]
        worst = max(selected_ranked, key=lambda item: (item[0], item[2].tautomer_id))
        selected.remove(worst[2].tautomer_id)
    selected.add(input_candidate.tautomer_id)
    return selected


def _single_candidate_rank(candidate: _Candidate, reason: str) -> list[_RankedCandidate]:
    return [
        _RankedCandidate(
            candidate=candidate,
            energy_kcal_mol=None,
            relative_energy_kcal_mol=None,
            auto3d_rank=1,
            selected=True,
            reason="selected_single_tautomer_candidate",
            source="rdkit_fallback",
            warnings=(reason,),
        )
    ]


def _fallback_rank(candidates: list[_Candidate], fallback_reason: str, config: RunConfig) -> list[_RankedCandidate]:
    selected = _fallback_selected_ids(candidates, config)
    scored = sorted(
        candidates,
        key=lambda candidate: (
            float("inf") if candidate.rdkit_score is None else -candidate.rdkit_score,
            not candidate.is_canonical_tautomer,
            not candidate.is_input_tautomer,
            candidate.isomeric_smiles,
        ),
    )
    rank_by_id = {candidate.tautomer_id: rank for rank, candidate in enumerate(scored, start=1)}
    return [
        _RankedCandidate(
            candidate=candidate,
            energy_kcal_mol=None,
            relative_energy_kcal_mol=None,
            auto3d_rank=rank_by_id.get(candidate.tautomer_id),
            selected=candidate.tautomer_id in selected,
            reason="selected_by_rdkit_fallback" if candidate.tautomer_id in selected else "rejected_by_rdkit_fallback",
            source="rdkit_fallback",
            warnings=(fallback_reason,),
        )
        for candidate in candidates
    ]


def _fallback_selected_ids(candidates: list[_Candidate], config: RunConfig) -> set[str]:
    top_n = max(1, min(config.tautomer_filtering.tauto_k, len(candidates)))
    selected: set[str] = set()
    for predicate in (
        lambda candidate: candidate.is_canonical_tautomer,
        lambda candidate: candidate.is_input_tautomer,
    ):
        candidate = next((item for item in candidates if predicate(item)), None)
        if candidate is not None:
            selected.add(candidate.tautomer_id)
    scored = sorted(
        candidates,
        key=lambda candidate: (
            float("inf") if candidate.rdkit_score is None else -candidate.rdkit_score,
            candidate.isomeric_smiles,
        ),
    )
    for candidate in scored:
        if len(selected) >= top_n:
            break
        selected.add(candidate.tautomer_id)
    if len(selected) > top_n:
        scored_selected = [candidate for candidate in scored if candidate.tautomer_id in selected]
        selected = {candidate.tautomer_id for candidate in scored_selected[:top_n]}
    return selected


def _records_from_ranked(protomer: ProtomerRecord, ranked: list[_RankedCandidate], config: RunConfig, output_dir: Path) -> list[TautomerRecord]:
    records: list[TautomerRecord] = []
    selected = sorted(
        ranked,
        key=lambda item: (
            float("inf") if item.relative_energy_kcal_mol is None else item.relative_energy_kcal_mol,
            item.candidate.isomeric_smiles,
        ),
    )
    for index, item in enumerate(selected, start=1):
        candidate = item.candidate
        metadata = {
            "auto3d_tautomer_filtering": {
                "selected": True,
                "source": item.source,
                "reason": item.reason,
                "rank": item.auto3d_rank,
                "energy_kcal_mol": item.energy_kcal_mol,
                "relative_energy_kcal_mol": item.relative_energy_kcal_mol,
                "score_is_population_estimate": False,
                "scope": "fast potential-energy tautomer filter before stereoisomer enumeration",
                "tauto_engine": config.tautomer_filtering.tauto_engine,
                "optimizing_engine": config.tautomer_filtering.optimizing_engine,
            }
        }
        records.append(
            TautomerRecord(
                id=candidate.tautomer_id,
                parent_id=protomer.id,
                input_molecule_id=protomer.input_molecule_id,
                molname=protomer.molname,
                canonical_smiles=candidate.canonical_smiles,
                isomeric_smiles=candidate.isomeric_smiles,
                molecular_formula=_formula(candidate.molecule),
                formal_charge=Chem.GetFormalCharge(candidate.molecule),
                explicit_proton_count=_explicit_proton_count(candidate.molecule),
                source_software=item.source,
                source_python_function="dsvr.chemistry.tautomer_auto3d_filter.filter_tautomers_with_auto3d",
                output_paths=[
                    output_dir / f"{protomer.id}_tautomers.sdf",
                    output_dir / f"{protomer.id}_tautomers.csv",
                    output_dir / "tautomers_selected.csv",
                    output_dir / "tautomers_rejected.csv",
                ],
                warnings=[
                    "Auto3D tautomer filtering is a fast potential-energy triage step, not a solution abundance estimate.",
                    *item.warnings,
                ],
                metadata=metadata,
                tautomer_index=index,
                rdkit_mol=candidate.molecule,
            )
        )
    return records


def _write_tautomer_sdf(path: Path, records: list[TautomerRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        molecule = Chem.Mol(record.rdkit_mol)
        molecule.SetProp("_Name", record.id)
        filtering = record.metadata.get("auto3d_tautomer_filtering", {})
        relative = filtering.get("relative_energy_kcal_mol")
        for key, value in {
            "DSVR_STAGE": record.stage_name,
            "DSVR_INPUT_ID": record.input_molecule_id,
            "DSVR_PARENT_PROTOMER_ID": record.parent_id or "",
            "DSVR_TAUTOMER_ID": record.id,
            "DSVR_MOLNAME": record.molname,
            "DSVR_CANONICAL_SMILES": record.canonical_smiles or "",
            "DSVR_ISOMERIC_SMILES": record.isomeric_smiles or "",
            "DSVR_FORMULA": record.molecular_formula or "",
            "DSVR_FORMAL_CHARGE": str(record.formal_charge),
            "DSVR_EXPLICIT_PROTON_COUNT": str(record.explicit_proton_count),
            "DSVR_TAUTOMER_FILTER_SOURCE": str(filtering.get("source", "")),
            "DSVR_TAUTOMER_FILTER_REASON": str(filtering.get("reason", "")),
            "E_tautomer_relative(kcal/mol)": "" if relative is None else str(relative),
        }.items():
            molecule.SetProp(key, value)
        writer.write(molecule)
    writer.close()


def _write_tautomer_csv(path: Path, records: list[TautomerRecord]) -> None:
    rows = []
    for record in records:
        filtering = record.metadata.get("auto3d_tautomer_filtering", {})
        rows.append(
            {
                "id": record.id,
                "parent_id": record.parent_id,
                "input_molecule_id": record.input_molecule_id,
                "molname": record.molname,
                "canonical_smiles": record.canonical_smiles,
                "isomeric_smiles": record.isomeric_smiles,
                "source": filtering.get("source"),
                "reason": filtering.get("reason"),
                "auto3d_rank": filtering.get("rank"),
                "energy_kcal_mol": filtering.get("energy_kcal_mol"),
                "relative_energy_kcal_mol": filtering.get("relative_energy_kcal_mol"),
                "warnings": " | ".join(record.warnings),
            }
        )
    _write_rows(path, rows)


def _candidate_row(protomer: ProtomerRecord, candidate: _Candidate, warning: str | None) -> dict[str, Any]:
    return {
        "protomer_id": protomer.id,
        "tautomer_id": candidate.tautomer_id,
        "input_molecule_id": protomer.input_molecule_id,
        "molname": protomer.molname,
        "canonical_smiles": candidate.canonical_smiles,
        "isomeric_smiles": candidate.isomeric_smiles,
        "rdkit_score": candidate.rdkit_score,
        "is_input_tautomer": candidate.is_input_tautomer,
        "is_canonical_tautomer": candidate.is_canonical_tautomer,
        "warning": warning or "",
    }


def _ranked_row(protomer: ProtomerRecord, ranked: _RankedCandidate) -> dict[str, Any]:
    candidate = ranked.candidate
    return {
        "protomer_id": protomer.id,
        "tautomer_id": candidate.tautomer_id,
        "input_molecule_id": protomer.input_molecule_id,
        "molname": protomer.molname,
        "canonical_smiles": candidate.canonical_smiles,
        "isomeric_smiles": candidate.isomeric_smiles,
        "selected": ranked.selected,
        "reason": ranked.reason,
        "source": ranked.source,
        "auto3d_rank": ranked.auto3d_rank,
        "energy_kcal_mol": ranked.energy_kcal_mol,
        "relative_energy_kcal_mol": ranked.relative_energy_kcal_mol,
        "rdkit_score": candidate.rdkit_score,
        "score_is_population_estimate": False,
        "warnings": " | ".join(ranked.warnings),
    }


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = [
        "protomer_id",
        "tautomer_id",
        "id",
        "parent_id",
        "input_molecule_id",
        "molname",
        "canonical_smiles",
        "isomeric_smiles",
        "selected",
        "reason",
        "source",
        "auto3d_rank",
        "energy_kcal_mol",
        "relative_energy_kcal_mol",
        "rdkit_score",
        "is_input_tautomer",
        "is_canonical_tautomer",
        "score_is_population_estimate",
        "warning",
        "warnings",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _score_tautomer(enumerator: rdMolStandardize.TautomerEnumerator, molecule: Chem.Mol) -> float | None:
    score = getattr(enumerator, "ScoreTautomer", None)
    if score is None:
        return None
    try:
        return float(score(molecule))
    except (RuntimeError, ValueError, TypeError):
        return None


def _formula(molecule: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(molecule)


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)
