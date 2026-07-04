from __future__ import annotations

import csv
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.config import RunConfig
from dsvr.models import (
    ProtomerRecord,
    SeedConformerRecord,
    StereoRecord,
    ThermoRecord,
    make_seed_id,
)
from dsvr.ranking.boltzmann import R_KCAL_MOL_K
from dsvr.runners.auto3d_runner import inspect_auto3d, run_auto3d
from dsvr.utils.units import hartree_to_kcal_mol


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
        mpi_np=config.seeding.auto3d_mpi_np,
        cpu_workers=config.seeding.auto3d_cpu_workers,
        memory_gb=config.seeding.auto3d_memory_gb,
        capacity=config.seeding.auto3d_capacity,
        max_confs=config.seeding.auto3d_max_confs,
        patience=config.seeding.auto3d_patience,
        threshold=config.seeding.auto3d_threshold,
        opt_steps=config.seeding.auto3d_opt_steps,
        use_gpu=config.seeding.auto3d_use_gpu,
        stream_output=config.logging.tail_subprocess_logs,
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


def generate_auto3d_seeds_from_protomers(
    protomer_records: list[ProtomerRecord],
    config: RunConfig,
) -> list[SeedConformerRecord]:
    output_dir = config.output_dir / "seeding" / "auto3d_protocol"
    output_dir.mkdir(parents=True, exist_ok=True)
    input_smi = output_dir / "auto3d_protomer_input.smi"
    _write_auto3d_protomer_input(input_smi, protomer_records)
    output_sdf, command = run_auto3d(
        input_smi,
        output_dir,
        k=config.seeding.auto3d_k,
        model=config.seeding.auto3d_model,
        internal_tautomer_stereo_enum=True,
        mpi_np=config.seeding.auto3d_mpi_np,
        cpu_workers=config.seeding.auto3d_cpu_workers,
        memory_gb=config.seeding.auto3d_memory_gb,
        capacity=config.seeding.auto3d_capacity,
        max_confs=config.seeding.auto3d_max_confs,
        patience=config.seeding.auto3d_patience,
        threshold=config.seeding.auto3d_threshold,
        opt_steps=config.seeding.auto3d_opt_steps,
        use_gpu=config.seeding.auto3d_use_gpu,
        stream_output=config.logging.tail_subprocess_logs,
    )
    records = _records_from_auto3d_protomer_output(
        output_sdf,
        protomer_records,
        config=config,
        command=command,
        output_dir=output_dir,
    )
    _write_auto3d_seed_sdf(output_dir / "auto3d_protocol_seeds.sdf", records)
    _write_auto3d_seed_csv(output_dir / "auto3d_protocol_seeds.csv", records)
    return records


def reduce_auto3d_entropy_ensemble(
    records: list[SeedConformerRecord],
    config: RunConfig,
) -> list[ThermoRecord]:
    output_dir = config.output_dir / "auto3d_entropy"
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[tuple[str, str, str, int | None, int | None], list[SeedConformerRecord]] = {}
    for record in records:
        if record.energy_kcal_mol is None:
            continue
        key = (
            record.input_molecule_id,
            record.isomeric_smiles or record.canonical_smiles or record.id,
            record.molecular_formula or "",
            record.formal_charge,
            record.explicit_proton_count,
        )
        grouped.setdefault(key, []).append(record)

    thermo_records = []
    for index, ensemble in enumerate(grouped.values(), start=1):
        sorted_ensemble = sorted(
            ensemble,
            key=lambda item: (item.energy_kcal_mol or float("inf"), item.id),
        )
        representative = sorted_ensemble[0]
        free_energy, entropy = _ensemble_free_energy_and_entropy(sorted_ensemble, config)
        metadata = {
            "auto3d_entropy": {
                "protocol": "molscrub_auto3d_entropy",
                "source_record_ids": [item.id for item in sorted_ensemble],
                "conformer_count": len(sorted_ensemble),
                "temperature_kelvin": config.chemistry.temperature_kelvin,
                "reduction": "configurational_partition_function",
                "energy_source": "Auto3D optimized conformer energies",
            }
        }
        thermo_records.append(
            ThermoRecord(
                id=f"{representative.input_molecule_id}_auto3d_entropy_{index:04d}",
                parent_id=representative.id,
                input_molecule_id=representative.input_molecule_id,
                molname=representative.molname,
                canonical_smiles=representative.canonical_smiles,
                isomeric_smiles=representative.isomeric_smiles,
                molecular_formula=representative.molecular_formula,
                formal_charge=representative.formal_charge,
                explicit_proton_count=representative.explicit_proton_count,
                source_software="auto3d",
                source_command=representative.source_command,
                source_python_function=(
                    "dsvr.chemistry.conformers_auto3d.reduce_auto3d_entropy_ensemble"
                ),
                output_paths=[
                    output_dir / "auto3d_entropy_records.csv",
                    output_dir / "auto3d_entropy_records.jsonl",
                ],
                warnings=[
                    "Auto3D entropy ranking uses a configurational conformer-ensemble "
                    "free energy from optimized Auto3D energies; it is not a full "
                    "vibrational/RRHO thermochemistry calculation."
                ],
                metadata=metadata,
                temperature_kelvin=config.chemistry.temperature_kelvin,
                free_energy_kcal_mol=free_energy,
                entropy_cal_mol_k=entropy,
            )
        )
    _write_auto3d_entropy_outputs(output_dir, thermo_records)
    return thermo_records


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


def _write_auto3d_protomer_input(path: Path, protomer_records: list[ProtomerRecord]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in protomer_records:
            smiles = record.isomeric_smiles or record.canonical_smiles or ""
            handle.write(f"{smiles} {record.id}\n")


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


def _records_from_auto3d_protomer_output(
    output_sdf: Path,
    protomer_records: list[ProtomerRecord],
    *,
    config: RunConfig,
    command: list[str],
    output_dir: Path,
) -> list[SeedConformerRecord]:
    protomer_by_id = {record.id: record for record in protomer_records}
    fallback_protomer = protomer_records[0] if len(protomer_records) == 1 else None
    supplier = Chem.SDMolSupplier(str(output_sdf), sanitize=True, removeHs=False)
    records: list[SeedConformerRecord] = []
    seen: set[tuple[str, str | None, str | None]] = set()
    for index, molecule in enumerate(supplier, start=1):
        if molecule is None:
            continue
        protomer_id = _source_protomer_id(molecule, protomer_by_id)
        parent = protomer_by_id.get(protomer_id or "") or fallback_protomer
        if parent is None:
            continue
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        energy, energy_prop, energy_units = _extract_energy_with_units(molecule)
        key = (
            parent.id,
            isomeric_smiles,
            None if energy is None else f"{energy:.8f}",
        )
        if key in seen:
            continue
        seen.add(key)
        metadata = {
            "auto3d": {
                "k": config.seeding.auto3d_k,
                "model": config.seeding.auto3d_model,
                "lineage_mode": "protomer_to_auto3d_internal_tautomer_stereo_enum",
                "internal_tautomer_stereo_enum": True,
                "source_protomer_id": protomer_id,
                "energy_property": energy_prop,
                "energy_units": energy_units,
                "max_confs": config.seeding.auto3d_max_confs,
                "mpi_np": config.seeding.auto3d_mpi_np,
                "cpu_workers": config.seeding.auto3d_cpu_workers,
                "memory_gb": config.seeding.auto3d_memory_gb,
                "capacity": config.seeding.auto3d_capacity,
                "patience": config.seeding.auto3d_patience,
                "threshold": config.seeding.auto3d_threshold,
                "opt_steps": config.seeding.auto3d_opt_steps,
                "use_gpu": config.seeding.auto3d_use_gpu,
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
                source_python_function=(
                    "dsvr.chemistry.conformers_auto3d.generate_auto3d_seeds_from_protomers"
                ),
                output_paths=[
                    output_sdf,
                    output_dir / "auto3d_protocol_seeds.sdf",
                    output_dir / "auto3d_protocol_seeds.csv",
                ],
                warnings=[
                    "Auto3D internally enumerated tautomers/stereoisomers after molscrub "
                    "protomer generation; lineage is anchored at the protomer."
                ],
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
    energy, _prop, _units = _extract_energy_with_units(molecule)
    return energy


def _extract_energy_with_units(molecule: Chem.Mol) -> tuple[float | None, str | None, str | None]:
    prop = _energy_property_name(molecule)
    if prop is None:
        return None, None, None
    try:
        value = float(molecule.GetProp(prop))
    except ValueError:
        return None, prop, None
    if prop == "E_tot":
        return hartree_to_kcal_mol(value), prop, "hartree_to_kcal_mol"
    return value, prop, "kcal_mol"


def _energy_property_name(molecule: Chem.Mol) -> str | None:
    for prop in ("E_kcal_mol", "E_tot", "E", "energy", "Energy", "E_rel", "E_relative"):
        if molecule.HasProp(prop):
            return prop
    return None


def _source_protomer_id(
    molecule: Chem.Mol,
    protomer_by_id: dict[str, ProtomerRecord],
) -> str | None:
    for prop in ("DSVR_PROTOMER_ID", "protomer_id", "parent_id", "_Name"):
        if not molecule.HasProp(prop):
            continue
        value = molecule.GetProp(prop).strip()
        if value in protomer_by_id:
            return value
        for protomer_id in protomer_by_id:
            if value.startswith(protomer_id):
                return protomer_id
    return None


def _ensemble_free_energy_and_entropy(
    records: list[SeedConformerRecord],
    config: RunConfig,
) -> tuple[float, float]:
    energies = [record.energy_kcal_mol for record in records if record.energy_kcal_mol is not None]
    if not energies:
        return 0.0, 0.0
    temperature = config.chemistry.temperature_kelvin
    rt = R_KCAL_MOL_K * temperature
    minimum = min(energies)
    weights = [math.exp(-(energy - minimum) / rt) for energy in energies]
    partition = sum(weights)
    free_energy = minimum - rt * math.log(partition)
    probabilities = [weight / partition for weight in weights]
    mean_energy = sum(
        probability * energy
        for probability, energy in zip(probabilities, energies, strict=True)
    )
    entropy_cal_mol_k = (mean_energy - free_energy) * 1000.0 / temperature
    return free_energy, entropy_cal_mol_k


def _write_auto3d_entropy_outputs(output_dir: Path, records: list[ThermoRecord]) -> None:
    csv_path = output_dir / "auto3d_entropy_records.csv"
    jsonl_path = output_dir / "auto3d_entropy_records.jsonl"
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
        "temperature_kelvin",
        "free_energy_kcal_mol",
        "entropy_cal_mol_k",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = record.model_dump(mode="json")
            writer.writerow({column: row.get(column) for column in columns})
    with jsonl_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(record.model_dump_json() + "\n")


def _formula(molecule: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(molecule)


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)
