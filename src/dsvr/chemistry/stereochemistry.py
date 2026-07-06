from __future__ import annotations

import csv
import multiprocessing as mp
import queue
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)

from dsvr.config import RunConfig
from dsvr.models import StereoRecord, TautomerRecord, make_stereo_id


STEREO_TIMEOUT_FALLBACK = "STEREO_TIMEOUT_FALLBACK"


def enumerate_stereoisomers(
    tautomer_record: TautomerRecord,
    config: RunConfig,
) -> list[StereoRecord]:
    max_isomers = _max_stereoisomers(config)
    options = StereoEnumerationOptions(
        tryEmbedding=(
            config.stereoisomer_filtering.try_embedding
            and config.enumeration.stereo_try_embedding
        ),
        onlyUnassigned=(
            config.stereoisomer_filtering.only_unassigned
            and config.enumeration.stereo_only_unassigned
        ),
        unique=config.enumeration.stereo_unique,
        maxIsomers=max_isomers,
        rand=config.enumeration.stereo_random_seed,
    )
    input_mol = Chem.Mol(tautomer_record.rdkit_mol)
    extra_warnings: list[str] = []
    try:
        raw_stereoisomers = _enumerate_with_timeout(
            input_mol,
            timeout_seconds=config.stereoisomer_filtering.timeout_seconds_per_tautomer,
            try_embedding=options.tryEmbedding,
            only_unassigned=options.onlyUnassigned,
            unique=options.unique,
            max_isomers=max_isomers,
            random_seed=options.rand,
        )
        if (
            options.tryEmbedding
            and len(raw_stereoisomers) < max_isomers
            and _has_potential_double_bond_stereo(input_mol)
        ):
            retry_stereoisomers = _enumerate_with_timeout(
                input_mol,
                timeout_seconds=config.stereoisomer_filtering.timeout_seconds_per_tautomer,
                try_embedding=False,
                only_unassigned=options.onlyUnassigned,
                unique=options.unique,
                max_isomers=max_isomers,
                random_seed=options.rand,
            )
            if len(retry_stereoisomers) > len(raw_stereoisomers):
                raw_stereoisomers = retry_stereoisomers
                extra_warnings.append(
                    "RDKit tryEmbedding under-enumerated potential double-bond stereo; "
                    "retried enumeration without embedding to preserve E/Z candidates."
                )
    except TimeoutError:
        raw_stereoisomers = [input_mol]
        extra_warnings.append(
            f"{STEREO_TIMEOUT_FALLBACK}: RDKit stereoisomer enumeration timeout; "
            "retained input stereo state"
        )
    except RuntimeError as exc:
        raw_stereoisomers = [input_mol]
        extra_warnings.append(f"RDKit stereoisomer enumeration failed; retained input state: {exc}")
    output_dir = config.output_dir / "enumeration" / "stereoisomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    records = _records_from_stereoisomers(
        tautomer_record,
        raw_stereoisomers,
        config=config,
        options=options,
        output_dir=output_dir,
        extra_warnings=extra_warnings,
    )
    _write_stereo_sdf(output_dir / f"{tautomer_record.id}_stereoisomers.sdf", records)
    _write_stereo_csv(output_dir / f"{tautomer_record.id}_stereoisomers.csv", records)
    return records


def _records_from_stereoisomers(
    tautomer_record: TautomerRecord,
    stereoisomers: list[Chem.Mol],
    *,
    config: RunConfig,
    options: StereoEnumerationOptions,
    output_dir: Path,
    extra_warnings: list[str] | None = None,
) -> list[StereoRecord]:
    seen: set[str] = set()
    unique_stereoisomers: list[Chem.Mol] = []
    for stereoisomer in stereoisomers:
        isomeric_smiles = Chem.MolToSmiles(stereoisomer, canonical=True, isomericSmiles=True)
        if isomeric_smiles in seen:
            continue
        seen.add(isomeric_smiles)
        unique_stereoisomers.append(stereoisomer)

    cap = _max_stereoisomers(config)
    hit_cap = len(unique_stereoisomers) >= cap and len(stereoisomers) >= cap
    limited_stereoisomers = unique_stereoisomers[:cap]
    records: list[StereoRecord] = []
    for index, stereoisomer in enumerate(limited_stereoisomers, start=1):
        canonical_smiles = Chem.MolToSmiles(stereoisomer, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(stereoisomer, canonical=True, isomericSmiles=True)
        formula = _formula(stereoisomer)
        charge = Chem.GetFormalCharge(stereoisomer)
        proton_count = _explicit_proton_count(stereoisomer)
        metadata = {
            "candidate_generation_only": True,
            "rdkit_stereo_parameters": _stereo_parameters(options),
            "stereochemical_smiles": isomeric_smiles,
            "dedupe_key": {"isomeric_smiles": isomeric_smiles},
        }
        warnings = [
            *(extra_warnings or []),
            "RDKit stereoisomer enumeration is candidate generation only; dominance "
            "ranking occurs later.",
            "tryEmbedding is a heuristic filter and can be computationally expensive.",
        ]
        if config.enumeration.stereo_only_unassigned:
            warnings.append("Assigned stereochemistry was preserved by default.")
        else:
            warnings.append(
                "All stereocenters were eligible for enumeration, including assigned centers."
            )
        if hit_cap:
            warnings.append(
                "stereoisomer candidate count reached max_stereoisomers_per_tautomer; "
                f"candidates were limited to {cap}"
            )
        records.append(
            StereoRecord(
                id=make_stereo_id(
                    tautomer_record.id,
                    index,
                    canonical_smiles,
                    isomeric_smiles,
                    metadata,
                ),
                parent_id=tautomer_record.id,
                input_molecule_id=tautomer_record.input_molecule_id,
                molname=tautomer_record.molname,
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                molecular_formula=formula,
                formal_charge=charge,
                explicit_proton_count=proton_count,
                source_software="rdkit",
                source_python_function="dsvr.chemistry.stereochemistry.enumerate_stereoisomers",
                output_paths=[
                    output_dir / f"{tautomer_record.id}_stereoisomers.sdf",
                    output_dir / f"{tautomer_record.id}_stereoisomers.csv",
                ],
                warnings=warnings,
                metadata=metadata,
                stereo_index=index,
                rdkit_mol=stereoisomer,
            )
        )
    return records


def _has_potential_double_bond_stereo(molecule: Chem.Mol) -> bool:
    for bond in molecule.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE or bond.IsInRing():
            continue
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetAtomicNum() == 6 and end.GetAtomicNum() == 6:
            if begin.GetDegree() > 1 and end.GetDegree() > 1:
                return True
    return False


def _max_stereoisomers(config: RunConfig) -> int:
    return min(
        config.stereoisomer_filtering.max_stereoisomers_per_tautomer,
        config.enumeration.max_stereoisomers_per_tautomer,
    )


def _enumerate_with_timeout(
    molecule: Chem.Mol,
    *,
    timeout_seconds: int,
    try_embedding: bool,
    only_unassigned: bool,
    unique: bool,
    max_isomers: int,
    random_seed: int,
) -> list[Chem.Mol]:
    output_queue: mp.Queue = mp.Queue(maxsize=1)
    process = mp.Process(
        target=_stereo_worker,
        args=(
            Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
            try_embedding,
            only_unassigned,
            unique,
            max_isomers,
            random_seed,
            output_queue,
        ),
    )
    try:
        process.start()
        process.join(timeout_seconds)
        if process.is_alive():
            process.terminate()
            process.join(2)
            if process.is_alive():
                process.kill()
                process.join(2)
            raise TimeoutError("RDKit stereoisomer enumeration timed out")
        try:
            payload = output_queue.get_nowait()
        except queue.Empty as exc:
            raise RuntimeError("RDKit stereoisomer worker produced no output") from exc
        if payload.get("status") != "ok":
            raise RuntimeError(str(payload.get("error", "unknown stereo worker error")))
        molecules = [
            mol
            for mol in (
                Chem.MolFromMolBlock(block, sanitize=True, removeHs=False)
                for block in payload.get("molblocks", [])
            )
            if mol is not None
        ]
        return molecules or [Chem.Mol(molecule)]
    finally:
        if process.is_alive():
            process.terminate()
            process.join(2)
            if process.is_alive():
                process.kill()
                process.join(2)
        output_queue.close()
        output_queue.join_thread()
        close = getattr(process, "close", None)
        if close is not None:
            close()


def _stereo_worker(
    smiles: str,
    try_embedding: bool,
    only_unassigned: bool,
    unique: bool,
    max_isomers: int,
    random_seed: int,
    output_queue: mp.Queue,
) -> None:
    try:
        molecule = Chem.MolFromSmiles(smiles, sanitize=True)
        if molecule is None:
            raise ValueError("could not parse stereo worker molecule")
        options = StereoEnumerationOptions(
            tryEmbedding=try_embedding,
            onlyUnassigned=only_unassigned,
            unique=unique,
            maxIsomers=max_isomers,
            rand=random_seed,
        )
        stereoisomers = list(EnumerateStereoisomers(molecule, options=options))
        output_queue.put(
            {
                "status": "ok",
                "molblocks": [
                    Chem.MolToMolBlock(stereoisomer)
                    for stereoisomer in stereoisomers[:max_isomers]
                ],
            }
        )
    except Exception as exc:  # pragma: no cover - exercised through parent process.
        output_queue.put({"status": "error", "error": str(exc)})


def read_tautomers_sdf(path: Path) -> list[TautomerRecord]:
    supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
    records: list[TautomerRecord] = []
    for index, molecule in enumerate(supplier, start=1):
        if molecule is None:
            continue
        tautomer_id = _prop_or_default(molecule, "DSVR_TAUTOMER_ID", f"tautomer_{index:06d}")
        input_id = _prop_or_default(molecule, "DSVR_INPUT_ID", tautomer_id)
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        records.append(
            TautomerRecord(
                id=tautomer_id,
                parent_id=_prop_or_default(molecule, "DSVR_PARENT_PROTOMER_ID", input_id),
                input_molecule_id=input_id,
                molname=_prop_or_default(molecule, "DSVR_MOLNAME", molecule.GetProp("_Name")),
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                molecular_formula=_formula(molecule),
                formal_charge=Chem.GetFormalCharge(molecule),
                explicit_proton_count=_explicit_proton_count(molecule),
                source_software="sdf",
                source_python_function="dsvr.chemistry.stereochemistry.read_tautomers_sdf",
                tautomer_index=index,
                rdkit_mol=molecule,
            )
        )
    return records


def _write_stereo_sdf(path: Path, records: list[StereoRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        mol = Chem.Mol(record.rdkit_mol)
        mol.SetProp("_Name", record.id)
        for key, value in {
            "DSVR_STAGE": record.stage_name,
            "DSVR_INPUT_ID": record.input_molecule_id,
            "DSVR_PARENT_TAUTOMER_ID": record.parent_id or "",
            "DSVR_STEREO_ID": record.id,
            "DSVR_MOLNAME": record.molname,
            "DSVR_CANONICAL_SMILES": record.canonical_smiles or "",
            "DSVR_ISOMERIC_SMILES": record.isomeric_smiles or "",
            "DSVR_STEREOCHEMICAL_SMILES": record.isomeric_smiles or "",
            "DSVR_FORMULA": record.molecular_formula or "",
            "DSVR_FORMAL_CHARGE": str(record.formal_charge),
            "DSVR_EXPLICIT_PROTON_COUNT": str(record.explicit_proton_count),
        }.items():
            mol.SetProp(key, value)
        writer.write(mol)
    writer.close()


def _write_stereo_csv(path: Path, records: list[StereoRecord]) -> None:
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


def _stereo_parameters(options: StereoEnumerationOptions) -> dict[str, int | bool]:
    return {
        "tryEmbedding": options.tryEmbedding,
        "onlyUnassigned": options.onlyUnassigned,
        "unique": options.unique,
        "maxIsomers": options.maxIsomers,
        "rand": options.rand,
    }


def _formula(molecule: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(molecule)


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)


def _prop_or_default(molecule: Chem.Mol, key: str, default: str) -> str:
    return molecule.GetProp(key) if molecule.HasProp(key) else default


def enumerate_stereoisomers_placeholder(
    smiles: str,
    max_stereoisomers: int = 64,
) -> list[str]:
    return [smiles][:max_stereoisomers]
