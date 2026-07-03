from __future__ import annotations

import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)

from dsvr.config import RunConfig
from dsvr.models import StereoRecord, TautomerRecord, make_stereo_id


def enumerate_stereoisomers(
    tautomer_record: TautomerRecord,
    config: RunConfig,
) -> list[StereoRecord]:
    options = StereoEnumerationOptions(
        tryEmbedding=config.enumeration.stereo_try_embedding,
        onlyUnassigned=config.enumeration.stereo_only_unassigned,
        unique=config.enumeration.stereo_unique,
        maxIsomers=config.enumeration.max_stereoisomers_per_tautomer,
        rand=config.enumeration.stereo_random_seed,
    )
    input_mol = Chem.Mol(tautomer_record.rdkit_mol)
    raw_stereoisomers = list(EnumerateStereoisomers(input_mol, options=options))
    output_dir = config.output_dir / "enumeration" / "stereoisomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    records = _records_from_stereoisomers(
        tautomer_record,
        raw_stereoisomers,
        config=config,
        options=options,
        output_dir=output_dir,
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
) -> list[StereoRecord]:
    seen: set[str] = set()
    unique_stereoisomers: list[Chem.Mol] = []
    for stereoisomer in stereoisomers:
        isomeric_smiles = Chem.MolToSmiles(stereoisomer, canonical=True, isomericSmiles=True)
        if isomeric_smiles in seen:
            continue
        seen.add(isomeric_smiles)
        unique_stereoisomers.append(stereoisomer)

    cap = config.enumeration.max_stereoisomers_per_tautomer
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
