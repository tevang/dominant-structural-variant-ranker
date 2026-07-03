from __future__ import annotations

import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.config import RunConfig
from dsvr.models import SeedConformerRecord, StereoRecord, make_seed_id
from dsvr.runners.auto3d_runner import inspect_auto3d, run_auto3d


def auto3d_available() -> bool:
    inspection = inspect_auto3d()
    return bool(
        inspection["python_api_available"]
        or inspection["auto3d_executable"]
        or inspection["auto3D_executable"]
        or inspection["Auto3D_executable"]
    )


def generate_auto3d_seeds(
    stereo_records: list[StereoRecord],
    config: RunConfig,
) -> list[SeedConformerRecord]:
    output_dir = config.output_dir / "seeding" / "auto3d"
    output_dir.mkdir(parents=True, exist_ok=True)
    input_sdf = output_dir / "auto3d_input.sdf"
    _write_auto3d_input(input_sdf, stereo_records)
    output_sdf, command = run_auto3d(
        input_sdf,
        output_dir,
        k=config.seeding.auto3d_k,
        model=config.seeding.auto3d_model,
        internal_tautomer_stereo_enum=config.seeding.auto3d_internal_tautomer_stereo_enum,
    )
    records = _records_from_auto3d_output(
        output_sdf,
        stereo_records,
        config=config,
        command=command,
        output_dir=output_dir,
    )
    _write_auto3d_seed_sdf(output_dir / "auto3d_seeds.sdf", records)
    _write_auto3d_seed_csv(output_dir / "auto3d_seeds.csv", records)
    return records


def _write_auto3d_input(path: Path, stereo_records: list[StereoRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in stereo_records:
        mol = Chem.Mol(record.rdkit_mol)
        mol.SetProp("_Name", record.id)
        mol.SetProp("DSVR_STEREO_ID", record.id)
        mol.SetProp("DSVR_INPUT_ID", record.input_molecule_id)
        mol.SetProp("DSVR_PARENT_TAUTOMER_ID", record.parent_id or "")
        mol.SetProp("DSVR_MOLNAME", record.molname)
        mol.SetProp("DSVR_CANONICAL_SMILES", record.canonical_smiles or "")
        mol.SetProp("DSVR_ISOMERIC_SMILES", record.isomeric_smiles or "")
        writer.write(mol)
    writer.close()


def _records_from_auto3d_output(
    output_sdf: Path,
    stereo_records: list[StereoRecord],
    *,
    config: RunConfig,
    command: list[str],
    output_dir: Path,
) -> list[SeedConformerRecord]:
    stereo_by_id = {record.id: record for record in stereo_records}
    fallback_stereo = stereo_records[0] if stereo_records else None
    supplier = Chem.SDMolSupplier(str(output_sdf), sanitize=True, removeHs=False)
    records: list[SeedConformerRecord] = []
    seen: set[tuple[str, str | None]] = set()
    for index, molecule in enumerate(supplier, start=1):
        if molecule is None:
            continue
        stereo_id = _source_stereo_id(molecule)
        parent = stereo_by_id.get(stereo_id or "")
        if parent is None and fallback_stereo is not None:
            parent = fallback_stereo
        if parent is None:
            continue
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        energy = _extract_energy(molecule)
        key = (isomeric_smiles, None if energy is None else f"{energy:.8f}")
        if key in seen:
            continue
        seen.add(key)
        lineage_mode = (
            "auto3d_internal_enum"
            if config.seeding.auto3d_internal_tautomer_stereo_enum
            else "post_stereo_seed"
        )
        warnings = []
        if config.seeding.auto3d_internal_tautomer_stereo_enum:
            warnings.append(
                "Auto3D internal tautomer/stereo enumeration was enabled; exact "
                "protomer-tautomer-stereo lineage may be less controlled."
            )
        else:
            warnings.append(
                "Auto3D internal tautomer/stereo enumeration disabled to avoid double enumeration."
            )
        metadata = {
            "auto3d": {
                "k": config.seeding.auto3d_k,
                "model": config.seeding.auto3d_model,
                "lineage_mode": lineage_mode,
                "internal_tautomer_stereo_enum": (
                    config.seeding.auto3d_internal_tautomer_stereo_enum
                ),
                "source_stereo_id": stereo_id,
                "energy_property": _energy_property_name(molecule),
            }
        }
        records.append(
            SeedConformerRecord(
                id=make_seed_id(parent.id, index, canonical_smiles, isomeric_smiles, metadata),
                parent_id=parent.id,
                input_molecule_id=parent.input_molecule_id,
                molname=parent.molname,
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                molecular_formula=_formula(molecule),
                formal_charge=Chem.GetFormalCharge(molecule),
                explicit_proton_count=_explicit_proton_count(molecule),
                source_software="auto3d",
                source_command=" ".join(command),
                source_python_function="dsvr.chemistry.conformers_auto3d.generate_auto3d_seeds",
                output_paths=[
                    output_sdf,
                    output_dir / "auto3d_seeds.sdf",
                    output_dir / "auto3d_seeds.csv",
                ],
                warnings=warnings,
                metadata=metadata,
                conformer_index=index,
                energy_kcal_mol=energy,
                rdkit_mol=molecule,
                rdkit_conformer_id=0 if molecule.GetNumConformers() else None,
                forcefield="auto3d",
                forcefield_status="auto3d_optimized" if energy is not None else "auto3d_output",
                minimization_converged=None,
                embedding_status="success",
            )
        )
    return records


def _write_auto3d_seed_sdf(path: Path, records: list[SeedConformerRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        mol = Chem.Mol(record.rdkit_mol)
        mol.SetProp("_Name", record.id)
        for key, value in {
            "DSVR_STAGE": record.stage_name,
            "DSVR_INPUT_ID": record.input_molecule_id,
            "DSVR_PARENT_STEREO_ID": record.parent_id or "",
            "DSVR_SEED_ID": record.id,
            "DSVR_MOLNAME": record.molname,
            "DSVR_CANONICAL_SMILES": record.canonical_smiles or "",
            "DSVR_ISOMERIC_SMILES": record.isomeric_smiles or "",
            "DSVR_FORMULA": record.molecular_formula or "",
            "DSVR_FORMAL_CHARGE": str(record.formal_charge),
            "DSVR_ENERGY_KCAL_MOL": (
                "" if record.energy_kcal_mol is None else str(record.energy_kcal_mol)
            ),
        }.items():
            mol.SetProp(key, value)
        writer.write(mol)
    writer.close()


def _write_auto3d_seed_csv(path: Path, records: list[SeedConformerRecord]) -> None:
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
        "conformer_index",
        "energy_kcal_mol",
        "forcefield_status",
        "embedding_status",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = record.model_dump(mode="json")
            row["warnings"] = " | ".join(record.warnings)
            writer.writerow({column: row.get(column) for column in columns})


def _source_stereo_id(molecule: Chem.Mol) -> str | None:
    for prop in ("DSVR_STEREO_ID", "stereo_id", "parent_id", "_Name"):
        if molecule.HasProp(prop):
            value = molecule.GetProp(prop).strip()
            if value:
                return value
    return None


def _extract_energy(molecule: Chem.Mol) -> float | None:
    prop = _energy_property_name(molecule)
    if prop is None:
        return None
    try:
        return float(molecule.GetProp(prop))
    except ValueError:
        return None


def _energy_property_name(molecule: Chem.Mol) -> str | None:
    for prop in ("E_tot", "E", "energy", "Energy", "E_kcal_mol"):
        if molecule.HasProp(prop):
            return prop
    return None


def _formula(molecule: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(molecule)


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)
