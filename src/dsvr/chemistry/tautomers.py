from __future__ import annotations

import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from dsvr.config import RunConfig
from dsvr.models import ProtomerRecord, TautomerRecord, make_tautomer_id


def enumerate_tautomers(
    protomer_record: ProtomerRecord,
    config: RunConfig,
) -> list[TautomerRecord]:
    enumerator = rdMolStandardize.TautomerEnumerator()
    max_tautomers = config.enumeration.max_tautomers_per_protomer
    enumerator.SetMaxTautomers(max_tautomers)

    input_mol = Chem.Mol(protomer_record.rdkit_mol)
    original_isomeric = Chem.MolToSmiles(input_mol, canonical=True, isomericSmiles=True)
    original_chiral_centers = _chiral_centers(input_mol)
    raw_tautomers = list(enumerator.Enumerate(input_mol))
    output_dir = config.output_dir / "enumeration" / "tautomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    records = _records_from_tautomers(
        protomer_record,
        raw_tautomers,
        config=config,
        enumerator=enumerator,
        original_isomeric=original_isomeric,
        original_chiral_centers=original_chiral_centers,
        output_dir=output_dir,
    )
    _write_tautomer_sdf(output_dir / f"{protomer_record.id}_tautomers.sdf", records)
    _write_tautomer_csv(output_dir / f"{protomer_record.id}_tautomers.csv", records)
    return records


def _records_from_tautomers(
    protomer_record: ProtomerRecord,
    tautomers: list[Chem.Mol],
    *,
    config: RunConfig,
    enumerator: rdMolStandardize.TautomerEnumerator,
    original_isomeric: str,
    original_chiral_centers: list[tuple[int, str]],
    output_dir: Path,
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

    cap = config.enumeration.max_tautomers_per_protomer
    hit_cap = len(unique_tautomers) >= cap and len(tautomers) >= cap
    limited_tautomers = unique_tautomers[:cap]
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
            "dedupe_key": {
                "canonical_smiles": canonical_smiles,
                "isomeric_smiles": isomeric_smiles,
            },
            "not_stability_ranking": True,
        }
        warnings = [
            "RDKit tautomer enumeration is candidate generation only; no tautomer "
            "stability ranking is implied."
        ]
        warnings.extend(_stereo_warnings(tautomer, original_isomeric, original_chiral_centers))
        if hit_cap:
            warnings.append(
                "tautomer candidate count reached max_tautomers_per_protomer; candidates "
                f"were limited to {cap}"
            )
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
