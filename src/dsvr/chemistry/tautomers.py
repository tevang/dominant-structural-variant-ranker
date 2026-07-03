from __future__ import annotations

import csv
import multiprocessing as mp
import queue
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from dsvr.config import RunConfig
from dsvr.filtering.variant_score import score_tautomer
from dsvr.models import ProtomerRecord, TautomerRecord, make_tautomer_id


class TautomerEnumerationTimeout(RuntimeError):
    """Raised when RDKit tautomer enumeration exceeds the configured wall time."""


@dataclass(frozen=True)
class _TautomerEnumerationResult:
    molblocks: list[str]
    hit_cap: bool
    elapsed_seconds: float
    worker_warnings: list[str]


@dataclass(frozen=True)
class _TautomerSettings:
    max_tautomers: int
    max_transforms: int
    timeout_seconds: int
    strategy: str
    remove_bond_stereo: bool
    remove_sp3_stereo: bool
    reassign_stereo: bool


@dataclass(frozen=True)
class _TautomerScoringRecord:
    rdkit_mol: Chem.Mol
    isomeric_smiles: str
    canonical_smiles: str
    formal_charge: int
    metadata: dict[str, Any]


def enumerate_tautomers(
    protomer_record: ProtomerRecord,
    config: RunConfig,
) -> list[TautomerRecord]:
    input_mol = Chem.Mol(protomer_record.rdkit_mol)
    original_isomeric = Chem.MolToSmiles(input_mol, canonical=True, isomericSmiles=True)
    original_chiral_centers = _chiral_centers(input_mol)
    settings = _settings_from_config(config)
    params = _cleanup_parameters(settings)
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    timeout_warning = None
    fallback = False
    try:
        result = _enumerate_molblocks_with_timeout(input_mol, settings)
        raw_tautomers = _mols_from_molblocks(result.molblocks)
        worker_warnings = result.worker_warnings
        hit_cap = result.hit_cap
        elapsed_seconds = result.elapsed_seconds
    except TautomerEnumerationTimeout:
        raw_tautomers = [input_mol]
        worker_warnings = ["tautomer enumeration timeout"]
        hit_cap = False
        elapsed_seconds = float(settings.timeout_seconds)
        timeout_warning = "tautomer enumeration timeout"
        fallback = True
    except RuntimeError as exc:
        raw_tautomers = [input_mol]
        worker_warnings = [f"tautomer enumeration failed; parent fallback used: {exc}"]
        hit_cap = False
        elapsed_seconds = 0.0
        fallback = True

    output_dir = config.output_dir / "enumeration" / "tautomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    records = _records_from_tautomers(
        protomer_record,
        raw_tautomers,
        config=config,
        enumerator=enumerator,
        settings=settings,
        original_isomeric=original_isomeric,
        original_chiral_centers=original_chiral_centers,
        output_dir=output_dir,
        hit_cap=hit_cap,
        worker_warnings=worker_warnings,
        elapsed_seconds=elapsed_seconds,
        fallback=fallback,
    )
    if timeout_warning is not None:
        records = [
            record.model_copy(update={"warnings": [*record.warnings, timeout_warning]})
            for record in records
        ]
    _write_tautomer_sdf(output_dir / f"{protomer_record.id}_tautomers.sdf", records)
    _write_tautomer_csv(output_dir / f"{protomer_record.id}_tautomers.csv", records)
    return records


def _records_from_tautomers(
    protomer_record: ProtomerRecord,
    tautomers: list[Chem.Mol],
    *,
    config: RunConfig,
    enumerator: rdMolStandardize.TautomerEnumerator,
    settings: _TautomerSettings,
    original_isomeric: str,
    original_chiral_centers: list[tuple[int, str]],
    output_dir: Path,
    hit_cap: bool,
    worker_warnings: list[str],
    elapsed_seconds: float,
    fallback: bool,
) -> list[TautomerRecord]:
    seen: set[tuple[str, str]] = set()
    unique_tautomers: list[Chem.Mol] = []
    for tautomer in tautomers:
        canonical_smiles = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=True)
        key = (canonical_smiles, isomeric_smiles)
        if key in seen:
            continue
        seen.add(key)
        unique_tautomers.append(tautomer)

    cap = settings.max_tautomers
    if len(unique_tautomers) > cap:
        hit_cap = True
    limited_tautomers = _select_tautomer_subset(unique_tautomers, cap, protomer_record)
    records: list[TautomerRecord] = []
    for index, tautomer in enumerate(limited_tautomers, start=1):
        canonical_smiles = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=True)
        formula = _formula(tautomer)
        charge = Chem.GetFormalCharge(tautomer)
        proton_count = _explicit_proton_count(tautomer)
        metadata = {
            "candidate_generation_only": True,
            "rdkit_tautomer_parameters": _tautomer_parameters(enumerator),
            "tautomer_strategy": settings.strategy,
            "elapsed_seconds": elapsed_seconds,
            "dedupe_key": {
                "canonical_smiles": canonical_smiles,
                "isomeric_smiles": isomeric_smiles,
            },
            "not_stability_ranking": True,
        }
        if fallback:
            metadata["fallback"] = True
            metadata["fallback_reason"] = "tautomer enumeration timeout or worker failure"
        warnings = [
            "RDKit tautomer enumeration is candidate generation only; no tautomer "
            "stability ranking is implied."
        ]
        if settings.strategy == "exhaustive":
            warnings.append(
                "tautomer_strategy=exhaustive can be very expensive and may generate "
                "large tautomer sets."
            )
        warnings.extend(worker_warnings)
        warnings.extend(_stereo_warnings(tautomer, original_isomeric, original_chiral_centers))
        if hit_cap:
            warnings.append(
                "tautomer candidate count reached max_tautomers_per_protomer; candidates "
                f"were limited to {cap} using SVPScore/RDKit tautomer heuristic priority"
            )
        if fallback:
            warnings.append("fallback tautomer candidate is the parent protomer itself")
        records.append(
            TautomerRecord(
                id=make_tautomer_id(
                    protomer_record.id,
                    index,
                    canonical_smiles,
                    isomeric_smiles,
                    metadata,
                ),
                parent_id=protomer_record.id,
                input_molecule_id=protomer_record.input_molecule_id,
                molname=protomer_record.molname,
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                molecular_formula=formula,
                formal_charge=charge,
                explicit_proton_count=proton_count,
                source_software="rdkit",
                source_python_function="dsvr.chemistry.tautomers.enumerate_tautomers",
                output_paths=[
                    output_dir / f"{protomer_record.id}_tautomers.sdf",
                    output_dir / f"{protomer_record.id}_tautomers.csv",
                ],
                warnings=warnings,
                metadata=metadata,
                tautomer_index=index,
                rdkit_mol=tautomer,
            )
        )
    return records


def _write_tautomer_sdf(path: Path, records: list[TautomerRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        mol = Chem.Mol(record.rdkit_mol)
        mol.SetProp("_Name", record.id)
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
            "DSVR_TAUTOMER_SCOPE": "candidate_generation_only_not_stability_ranking",
            "DSVR_TAUTOMER_FALLBACK": str(bool(record.metadata.get("fallback", False))),
        }.items():
            mol.SetProp(key, value)
        writer.write(mol)
    writer.close()


def _write_tautomer_csv(path: Path, records: list[TautomerRecord]) -> None:
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
        "source_python_function",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = record.model_dump(mode="json")
            row["warnings"] = " | ".join(record.warnings)
            writer.writerow({column: row.get(column) for column in columns})


def read_protomers_sdf(path: Path) -> list[ProtomerRecord]:
    supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
    records: list[ProtomerRecord] = []
    for index, molecule in enumerate(supplier, start=1):
        if molecule is None:
            continue
        protomer_id = _prop_or_default(molecule, "DSVR_PROTOMER_ID", f"protomer_{index:06d}")
        input_id = _prop_or_default(molecule, "DSVR_INPUT_ID", protomer_id)
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        records.append(
            ProtomerRecord(
                id=protomer_id,
                parent_id=_prop_or_default(molecule, "DSVR_PARENT_ID", input_id),
                input_molecule_id=input_id,
                molname=_prop_or_default(molecule, "DSVR_MOLNAME", molecule.GetProp("_Name")),
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                molecular_formula=_formula(molecule),
                formal_charge=Chem.GetFormalCharge(molecule),
                explicit_proton_count=_explicit_proton_count(molecule),
                source_software="sdf",
                source_python_function="dsvr.chemistry.tautomers.read_protomers_sdf",
                protomer_index=index,
                rdkit_mol=molecule,
            )
        )
    return records


def _tautomer_parameters(enumerator: rdMolStandardize.TautomerEnumerator) -> dict[str, int | bool]:
    return {
        "max_tautomers": enumerator.GetMaxTautomers(),
        "max_transforms": enumerator.GetMaxTransforms(),
        "remove_sp3_stereo": enumerator.GetRemoveSp3Stereo(),
        "remove_bond_stereo": enumerator.GetRemoveBondStereo(),
        "reassign_stereo": enumerator.GetReassignStereo(),
    }


def _settings_from_config(config: RunConfig) -> _TautomerSettings:
    return _TautomerSettings(
        max_tautomers=config.enumeration.max_tautomers_per_protomer,
        max_transforms=config.enumeration.max_tautomer_transforms,
        timeout_seconds=config.enumeration.tautomer_timeout_seconds,
        strategy=config.enumeration.tautomer_strategy,
        remove_bond_stereo=config.enumeration.tautomer_remove_bond_stereo,
        remove_sp3_stereo=config.enumeration.tautomer_remove_sp3_stereo,
        reassign_stereo=config.enumeration.tautomer_reassign_stereo,
    )


def _cleanup_parameters(settings: _TautomerSettings) -> rdMolStandardize.CleanupParameters:
    params = rdMolStandardize.CleanupParameters()
    params.maxTautomers = settings.max_tautomers
    params.maxTransforms = settings.max_transforms
    params.tautomerRemoveBondStereo = settings.remove_bond_stereo
    params.tautomerRemoveSp3Stereo = settings.remove_sp3_stereo
    params.tautomerReassignStereo = settings.reassign_stereo
    return params


def _enumerate_molblocks_with_timeout(
    molecule: Chem.Mol,
    settings: _TautomerSettings,
) -> _TautomerEnumerationResult:
    output_queue: Any = mp.Queue()
    molblock = Chem.MolToMolBlock(molecule)
    started = time.monotonic()
    process = mp.Process(
        target=_tautomer_worker,
        args=(molblock, settings, output_queue),
    )
    process.start()
    process.join(settings.timeout_seconds)
    elapsed_seconds = time.monotonic() - started
    if process.is_alive():
        process.terminate()
        process.join(2)
        if process.is_alive():
            process.kill()
            process.join(2)
        raise TautomerEnumerationTimeout("tautomer enumeration timeout")
    try:
        payload = output_queue.get_nowait()
    except queue.Empty as exc:
        raise RuntimeError("tautomer worker produced no result") from exc
    if payload.get("status") != "ok":
        raise RuntimeError(str(payload.get("error", "unknown tautomer worker error")))
    return _TautomerEnumerationResult(
        molblocks=list(payload["molblocks"]),
        hit_cap=bool(payload["hit_cap"]),
        elapsed_seconds=elapsed_seconds,
        worker_warnings=list(payload.get("warnings", [])),
    )


def _tautomer_worker(molblock: str, settings: _TautomerSettings, output_queue: Any) -> None:
    try:
        molecule = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
        if molecule is None:
            raise ValueError("could not deserialize protomer MolBlock")
        enumerator = rdMolStandardize.TautomerEnumerator(_cleanup_parameters(settings))
        tautomers = list(enumerator.Enumerate(molecule))
        output_queue.put(
            {
                "status": "ok",
                "molblocks": [Chem.MolToMolBlock(tautomer) for tautomer in tautomers],
                "hit_cap": len(tautomers) >= settings.max_tautomers,
                "warnings": [],
            }
        )
    except Exception as exc:  # pragma: no cover - exercised through parent error path.
        output_queue.put({"status": "error", "error": str(exc)})


def _mols_from_molblocks(molblocks: list[str]) -> list[Chem.Mol]:
    molecules = []
    for molblock in molblocks:
        molecule = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
        if molecule is not None:
            molecules.append(molecule)
    return molecules


def _select_tautomer_subset(
    tautomers: list[Chem.Mol],
    cap: int,
    protomer_record: ProtomerRecord,
) -> list[Chem.Mol]:
    if len(tautomers) <= cap:
        return tautomers
    parent_aromatic = sum(
        1 for atom in protomer_record.rdkit_mol.GetAtoms() if atom.GetIsAromatic()
    )
    scored = []
    for index, tautomer in enumerate(tautomers):
        canonical = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=False)
        isomeric = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=True)
        record = _TautomerScoringRecord(
            rdkit_mol=tautomer,
            canonical_smiles=canonical,
            isomeric_smiles=isomeric,
            formal_charge=Chem.GetFormalCharge(tautomer),
            metadata={"parent_aromatic_atom_count": parent_aromatic},
        )
        penalty, _ = score_tautomer(record)
        scored.append((penalty, isomeric, index, tautomer))
    return [tautomer for _, _, _, tautomer in sorted(scored)[:cap]]


def _stereo_warnings(
    tautomer: Chem.Mol,
    original_isomeric: str,
    original_chiral_centers: list[tuple[int, str]],
) -> list[str]:
    warnings: list[str] = []
    tautomer_isomeric = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=True)
    tautomer_nonisomeric = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=False)
    original_nonisomeric = Chem.MolToSmiles(
        Chem.MolFromSmiles(original_isomeric),
        canonical=True,
        isomericSmiles=False,
    )
    if original_isomeric != tautomer_isomeric and original_nonisomeric == tautomer_nonisomeric:
        warnings.append(
            "RDKit tautomer enumeration changed isomeric SMILES without changing "
            "non-isomeric connectivity; stereo labels may have changed."
        )
    tautomer_chiral_centers = _chiral_centers(tautomer)
    if original_chiral_centers and tautomer_chiral_centers != original_chiral_centers:
        warnings.append(
            "RDKit tautomer enumeration changed assigned chiral centers; stereoisomer "
            "enumeration must occur after tautomer enumeration."
        )
    return warnings


def _chiral_centers(molecule: Chem.Mol) -> list[tuple[int, str]]:
    return Chem.FindMolChiralCenters(
        molecule,
        includeUnassigned=True,
        useLegacyImplementation=False,
    )


def _formula(molecule: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(molecule)


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)


def _prop_or_default(molecule: Chem.Mol, key: str, default: str) -> str:
    return molecule.GetProp(key) if molecule.HasProp(key) else default


def enumerate_tautomers_placeholder(smiles: str, max_tautomers: int = 64) -> list[str]:
    return [smiles][:max_tautomers]
