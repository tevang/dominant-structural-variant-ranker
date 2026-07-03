from __future__ import annotations

import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.config import RunConfig
from dsvr.models import (
    MoleculeInput,
    ProtomerRecord,
    make_protomer_id,
)
from dsvr.runners.molscrub_runner import generate_molscrub_candidates


def describe_protonation_scope(ph: float) -> str:
    return f"candidate-generation pH={ph:g}; no rigorous pH population correction"


def generate_protomer_candidates(
    mol_record: MoleculeInput,
    config: RunConfig,
) -> list[ProtomerRecord]:
    ph_low = (
        config.chemistry.ph_low
        if config.chemistry.ph_low is not None
        else config.chemistry.ph
    )
    ph_high = (
        config.chemistry.ph_high
        if config.chemistry.ph_high is not None
        else config.chemistry.ph
    )
    raw_candidates, source_software, source_command = generate_molscrub_candidates(
        mol_record.rdkit_mol,
        ph_low=ph_low,
        ph_high=ph_high,
    )
    output_dir = config.output_dir / "enumeration" / "protomers"
    output_dir.mkdir(parents=True, exist_ok=True)
    records = _records_from_candidates(
        mol_record,
        raw_candidates,
        config=config,
        source_software=source_software,
        source_command=source_command,
        output_dir=output_dir,
    )
    _write_protomer_sdf(output_dir / f"{mol_record.input_id}_protomers.sdf", records)
    _write_protomer_csv(output_dir / f"{mol_record.input_id}_protomers.csv", records)
    return records


def _records_from_candidates(
    mol_record: MoleculeInput,
    candidates: list[Chem.Mol],
    *,
    config: RunConfig,
    source_software: str,
    source_command: str,
    output_dir: Path,
) -> list[ProtomerRecord]:
    seen: set[tuple[str, int, str]] = set()
    unique_candidates: list[Chem.Mol] = []
    for candidate in candidates:
        canonical_smiles = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=False)
        key = (
            _formula(candidate),
            Chem.GetFormalCharge(candidate),
            canonical_smiles,
        )
        if key in seen:
            continue
        seen.add(key)
        unique_candidates.append(candidate)

    cap = config.enumeration.max_protomers_per_molecule
    hit_cap = len(unique_candidates) > cap
    limited_candidates = unique_candidates[:cap]
    records: list[ProtomerRecord] = []
    for index, candidate in enumerate(limited_candidates, start=1):
        canonical_smiles = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True)
        formula = _formula(candidate)
        charge = Chem.GetFormalCharge(candidate)
        proton_count = _explicit_proton_count(candidate)
        metadata = {
            "ph_low": config.chemistry.ph_low,
            "ph_high": config.chemistry.ph_high,
            "solvent": config.chemistry.solvent,
            "candidate_generation_only": True,
            "dedupe_key": {
                "formula": formula,
                "formal_charge": charge,
                "canonical_smiles": canonical_smiles,
            },
        }
        warnings = [
            "molscrub pH influence is candidate generation only; no rigorous pH population "
            "prediction is implied."
        ]
        if hit_cap:
            warnings.append(
                "protomer candidate count exceeded max_protomers_per_molecule; candidates "
                f"were truncated to {cap}"
            )
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
                ],
                warnings=warnings,
                metadata=metadata,
                protomer_index=index,
                rdkit_mol=candidate,
            )
        )
    return records


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
            "DSVR_PH_SCOPE": "candidate_generation_only",
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


def _formula(molecule: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(molecule)


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)
