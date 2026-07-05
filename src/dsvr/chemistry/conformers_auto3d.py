from __future__ import annotations

import csv
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.chemistry.conformers_rdkit import generate_rdkit_seeds
from dsvr.config import RunConfig
from dsvr.filtering.variant_score import score_variant
from dsvr.models import (
    ProtomerRecord,
    SeedConformerRecord,
    StereoRecord,
    ThermoRecord,
    make_seed_id,
)
from dsvr.ranking.boltzmann import R_KCAL_MOL_K
from dsvr.runners.auto3d_runner import Auto3DExecutionError, inspect_auto3d, run_auto3d
from dsvr.utils.units import hartree_to_kcal_mol


@dataclass(frozen=True)
class _Auto3DProtomerFeatures:
    protomer_id: str
    input_molecule_id: str
    molname: str
    heavy_atoms: int
    rotatable_bonds: int
    stereo_centers: int
    hetero_atoms: int
    formal_charge: int
    estimated_enum_pressure: int


@dataclass(frozen=True)
class _Auto3DBatchSettings:
    internal_tautomer_stereo_enum: bool
    mpi_np: int | None
    cpu_workers: int | None
    memory_gb: int | None
    capacity: int | None
    max_confs: int | None
    patience: int | None
    threshold: float | None
    opt_steps: int | None
    use_gpu: bool
    reason: str


@dataclass(frozen=True)
class _Auto3DProtomerBatch:
    index: int
    records: list[ProtomerRecord]
    settings: _Auto3DBatchSettings
    large_molecule: bool
    features: list[_Auto3DProtomerFeatures]


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

    batches = _plan_auto3d_protomer_batches(protomer_records, config)
    _write_auto3d_adaptive_plan(output_dir / "auto3d_adaptive_plan.csv", batches)

    records: list[SeedConformerRecord] = []
    for batch in batches:
        batch_dir = output_dir / f"batch_{batch.index:03d}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        batch_input_smi = batch_dir / f"auto3d_protomer_input_batch_{batch.index:03d}.smi"
        _write_auto3d_protomer_input(batch_input_smi, batch.records)
        settings = batch.settings
        output_sdf, command = run_auto3d(
            batch_input_smi,
            batch_dir,
            k=config.seeding.auto3d_k,
            model=config.seeding.auto3d_model,
            internal_tautomer_stereo_enum=settings.internal_tautomer_stereo_enum,
            mpi_np=settings.mpi_np,
            cpu_workers=settings.cpu_workers,
            memory_gb=settings.memory_gb,
            capacity=settings.capacity,
            max_confs=settings.max_confs,
            patience=settings.patience,
            threshold=settings.threshold,
            opt_steps=settings.opt_steps,
            use_gpu=settings.use_gpu,
            stream_output=config.logging.tail_subprocess_logs,
        )
        batch_records = _records_from_auto3d_protomer_output(
            output_sdf,
            batch.records,
            config=config,
            command=command,
            output_dir=output_dir,
            settings=settings,
        )
        batch_records = _fill_missing_auto3d_protomer_outputs_with_rdkit(
            batch.records,
            batch_records,
            config=config,
            command=command,
            output_sdf=output_sdf,
            output_dir=batch_dir,
            final_output_dir=output_dir,
        )
        records.extend(batch_records)

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


def score_auto3d_representative_variants(
    records: list[SeedConformerRecord],
    config: RunConfig,
) -> list[ThermoRecord]:
    output_dir = config.output_dir / "auto3d_representatives"
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[tuple[str, str, str, int | None, int | None], list[SeedConformerRecord]] = {}
    for record in records:
        grouped.setdefault(_auto3d_variant_key(record), []).append(record)

    representatives = [
        sorted(
            ensemble,
            key=lambda item: (
                float("inf") if item.energy_kcal_mol is None else item.energy_kcal_mol,
                item.id,
            ),
        )[0]
        for ensemble in grouped.values()
        if ensemble
    ]
    best_energy_by_population_group: dict[tuple[object, ...], float] = {}
    for representative in representatives:
        if representative.energy_kcal_mol is None:
            continue
        key = _plausibility_population_group_key(representative, config)
        best_energy_by_population_group[key] = min(
            best_energy_by_population_group.get(key, float("inf")),
            representative.energy_kcal_mol,
        )

    scored_records = []
    for index, representative in enumerate(representatives, start=1):
        ensemble = grouped[_auto3d_variant_key(representative)]
        best_energy = best_energy_by_population_group.get(
            _plausibility_population_group_key(representative, config)
        )
        breakdown = score_variant(
            representative,
            {
                "ph": config.chemistry.ph,
                "best_energy": best_energy,
                "achiral_environment": not config.stereo_filtering.solvent_is_chiral,
            },
        )
        breakdown_dict = asdict(breakdown)
        metadata = {
            "auto3d_representative": {
                "protocol": "molscrub_auto3d_representative",
                "source_record_ids": [item.id for item in ensemble],
                "representative_record_id": representative.id,
                "candidate_conformer_count": len(ensemble),
                "representative_rule": "lowest_auto3d_energy_then_stable_id",
                "auto3d_energy_kcal_mol": representative.energy_kcal_mol,
                "score_field": "free_energy_kcal_mol stores SVPScore plausibility penalty",
            },
            "svp_score": breakdown_dict,
        }
        warnings = sorted(
            set(
                [
                    *representative.warnings,
                    *breakdown.warnings,
                    "Auto3D representative protocol ranks variants by SVPScore plausibility; "
                    "this is not conformer entropy, RRHO thermochemistry, or QM rescoring.",
                    "The score is stored in free_energy_kcal_mol only to reuse the existing "
                    "ranking/population output schema.",
                ]
            )
        )
        scored_records.append(
            ThermoRecord(
                id=f"{representative.input_molecule_id}_auto3d_repr_{index:04d}",
                parent_id=representative.id,
                input_molecule_id=representative.input_molecule_id,
                molname=representative.molname,
                canonical_smiles=representative.canonical_smiles,
                isomeric_smiles=representative.isomeric_smiles,
                molecular_formula=representative.molecular_formula,
                formal_charge=representative.formal_charge,
                explicit_proton_count=representative.explicit_proton_count,
                source_software="dsvr-svpscore",
                source_command=representative.source_command,
                source_python_function=(
                    "dsvr.chemistry.conformers_auto3d."
                    "score_auto3d_representative_variants"
                ),
                output_paths=[
                    output_dir / "auto3d_representative_scores.csv",
                    output_dir / "auto3d_representative_scores.jsonl",
                ],
                warnings=warnings,
                metadata=metadata,
                temperature_kelvin=config.chemistry.temperature_kelvin,
                free_energy_kcal_mol=breakdown.total,
                entropy_cal_mol_k=None,
            )
        )
    _write_auto3d_representative_outputs(output_dir, scored_records)
    return scored_records


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


def _plan_auto3d_protomer_batches(
    protomer_records: list[ProtomerRecord],
    config: RunConfig,
) -> list[_Auto3DProtomerBatch]:
    if not protomer_records:
        return []

    small_records: list[ProtomerRecord] = []
    small_features: list[_Auto3DProtomerFeatures] = []
    large_batches: list[tuple[ProtomerRecord, _Auto3DProtomerFeatures]] = []
    for record in protomer_records:
        features = _auto3d_protomer_features(record)
        if _is_large_auto3d_protomer(features):
            large_batches.append((record, features))
        else:
            small_records.append(record)
            small_features.append(features)

    batches: list[_Auto3DProtomerBatch] = []
    batch_index = 1
    if small_records:
        batches.append(
            _Auto3DProtomerBatch(
                index=batch_index,
                records=small_records,
                settings=_default_auto3d_protocol_settings(config),
                large_molecule=False,
                features=small_features,
            )
        )
        batch_index += 1

    for record, features in large_batches:
        batches.append(
            _Auto3DProtomerBatch(
                index=batch_index,
                records=[record],
                settings=_large_auto3d_protocol_settings(features, config),
                large_molecule=True,
                features=[features],
            )
        )
        batch_index += 1
    return batches


def _auto3d_protomer_features(record: ProtomerRecord) -> _Auto3DProtomerFeatures:
    molecule = Chem.Mol(record.rdkit_mol)
    heavy_atoms = molecule.GetNumHeavyAtoms()
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(molecule)
    stereo_centers = len(Chem.FindMolChiralCenters(molecule, includeUnassigned=True))
    hetero_atoms = sum(1 for atom in molecule.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
    formal_charge = Chem.GetFormalCharge(molecule)
    estimated_enum_pressure = (
        max(1, hetero_atoms)
        * max(1, stereo_centers + 1)
        * max(1, rotatable_bonds // 2)
    )
    return _Auto3DProtomerFeatures(
        protomer_id=record.id,
        input_molecule_id=record.input_molecule_id,
        molname=record.molname,
        heavy_atoms=heavy_atoms,
        rotatable_bonds=rotatable_bonds,
        stereo_centers=stereo_centers,
        hetero_atoms=hetero_atoms,
        formal_charge=formal_charge,
        estimated_enum_pressure=estimated_enum_pressure,
    )


def _is_large_auto3d_protomer(features: _Auto3DProtomerFeatures) -> bool:
    return (
        features.heavy_atoms >= 45
        or features.rotatable_bonds >= 12
        or features.estimated_enum_pressure >= 120
    )


def _default_auto3d_protocol_settings(config: RunConfig) -> _Auto3DBatchSettings:
    return _Auto3DBatchSettings(
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
        reason="configured Auto3D protocol settings",
    )


def _large_auto3d_protocol_settings(
    features: _Auto3DProtomerFeatures,
    config: RunConfig,
) -> _Auto3DBatchSettings:
    available_cpus = os.cpu_count() or config.seeding.auto3d_mpi_np
    configured_mpi = config.seeding.auto3d_mpi_np or available_cpus
    active_workers = max(
        2,
        min(configured_mpi, available_cpus, max(4, features.rotatable_bonds // 2)),
    )
    target_threads_per_worker = 4 if config.seeding.auto3d_use_gpu else 2
    cpu_workers = max(2, math.ceil(active_workers / target_threads_per_worker))
    cpu_workers = min(cpu_workers, active_workers)
    if config.seeding.auto3d_cpu_workers is not None:
        cpu_workers = max(1, min(cpu_workers, config.seeding.auto3d_cpu_workers))

    max_confs = config.seeding.auto3d_max_confs
    if max_confs is None or max_confs > 1:
        max_confs = 1
    capacity = config.seeding.auto3d_capacity
    if capacity is None or capacity > max_confs:
        capacity = max_confs
    opt_steps = config.seeding.auto3d_opt_steps
    opt_steps = 250 if opt_steps is None else min(opt_steps, 250)
    patience = config.seeding.auto3d_patience
    patience = 30 if patience is None else min(patience, 30)
    memory_gb = config.seeding.auto3d_memory_gb
    memory_gb = 2 if memory_gb is None else max(memory_gb, 2)

    reason = (
        "large/high-enumeration protomer: disabled Auto3D internal tautomer/stereo "
        "enumeration, capped conformers, and sized MPI/thread groups to expected "
        "active conformer work rather than total CPU count"
    )
    return _Auto3DBatchSettings(
        internal_tautomer_stereo_enum=False,
        mpi_np=active_workers,
        cpu_workers=cpu_workers,
        memory_gb=memory_gb,
        capacity=capacity,
        max_confs=max_confs,
        patience=patience,
        threshold=config.seeding.auto3d_threshold,
        opt_steps=opt_steps,
        use_gpu=config.seeding.auto3d_use_gpu,
        reason=reason,
    )


def _write_auto3d_adaptive_plan(path: Path, batches: list[_Auto3DProtomerBatch]) -> None:
    fieldnames = [
        "batch_index",
        "protomer_id",
        "input_molecule_id",
        "molname",
        "large_molecule",
        "heavy_atoms",
        "rotatable_bonds",
        "stereo_centers",
        "hetero_atoms",
        "formal_charge",
        "estimated_enum_pressure",
        "internal_tautomer_stereo_enum",
        "mpi_np",
        "cpu_workers",
        "memory_gb",
        "capacity",
        "max_confs",
        "patience",
        "threshold",
        "opt_steps",
        "use_gpu",
        "reason",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for batch in batches:
            settings = batch.settings
            for features in batch.features:
                writer.writerow(
                    {
                        "batch_index": batch.index,
                        "protomer_id": features.protomer_id,
                        "input_molecule_id": features.input_molecule_id,
                        "molname": features.molname,
                        "large_molecule": batch.large_molecule,
                        "heavy_atoms": features.heavy_atoms,
                        "rotatable_bonds": features.rotatable_bonds,
                        "stereo_centers": features.stereo_centers,
                        "hetero_atoms": features.hetero_atoms,
                        "formal_charge": features.formal_charge,
                        "estimated_enum_pressure": features.estimated_enum_pressure,
                        "internal_tautomer_stereo_enum": settings.internal_tautomer_stereo_enum,
                        "mpi_np": settings.mpi_np,
                        "cpu_workers": settings.cpu_workers,
                        "memory_gb": settings.memory_gb,
                        "capacity": settings.capacity,
                        "max_confs": settings.max_confs,
                        "patience": settings.patience,
                        "threshold": settings.threshold,
                        "opt_steps": settings.opt_steps,
                        "use_gpu": settings.use_gpu,
                        "reason": settings.reason,
                    }
                )


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
    settings: _Auto3DBatchSettings | None = None,
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
        batch_settings = settings or _default_auto3d_protocol_settings(config)
        lineage_mode = (
            "protomer_to_auto3d_internal_tautomer_stereo_enum"
            if batch_settings.internal_tautomer_stereo_enum
            else "protomer_direct_auto3d_large_molecule_constrained"
        )
        metadata = {
            "auto3d": {
                "k": config.seeding.auto3d_k,
                "model": config.seeding.auto3d_model,
                "lineage_mode": lineage_mode,
                "internal_tautomer_stereo_enum": batch_settings.internal_tautomer_stereo_enum,
                "source_protomer_id": protomer_id,
                "energy_property": energy_prop,
                "energy_units": energy_units,
                "max_confs": batch_settings.max_confs,
                "mpi_np": batch_settings.mpi_np,
                "cpu_workers": batch_settings.cpu_workers,
                "memory_gb": batch_settings.memory_gb,
                "capacity": batch_settings.capacity,
                "patience": batch_settings.patience,
                "threshold": batch_settings.threshold,
                "opt_steps": batch_settings.opt_steps,
                "use_gpu": batch_settings.use_gpu,
                "adaptive_reason": batch_settings.reason,
            }
        }
        if batch_settings.internal_tautomer_stereo_enum:
            warnings = [
                "Auto3D internally enumerated tautomers/stereoisomers after molscrub "
                "protomer generation; lineage is anchored at the protomer."
            ]
        else:
            warnings = [
                "Auto3D internal tautomer/stereo enumeration was disabled for this "
                "large/high-enumeration protomer to avoid excessive conformer expansion; "
                "lineage remains anchored at the molscrub protomer."
            ]
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


def _fill_missing_auto3d_protomer_outputs_with_rdkit(
    protomer_records: list[ProtomerRecord],
    records: list[SeedConformerRecord],
    *,
    config: RunConfig,
    command: list[str],
    output_sdf: Path,
    output_dir: Path,
    final_output_dir: Path | None = None,
) -> list[SeedConformerRecord]:
    protomer_by_id = {record.id: record for record in protomer_records}
    observed_ids = {record.parent_id for record in records if record.parent_id}
    missing_ids = sorted(protomer_by_id.keys() - observed_ids)
    if not missing_ids:
        return records

    fallback_records: list[SeedConformerRecord] = []
    fallback_failures: list[str] = []
    fallback_reason = _auto3d_partial_failure_reason(output_sdf, output_dir)
    for protomer_id in missing_ids:
        protomer = protomer_by_id[protomer_id]
        fallback = _rdkit_fallback_seed_for_missing_protomer(
            protomer,
            config=config,
            command=command,
            output_sdf=output_sdf,
            output_dir=final_output_dir or output_dir,
            reason=fallback_reason,
        )
        if fallback is None:
            fallback_failures.append(protomer_id)
            continue
        fallback_records.append(fallback)

    if fallback_failures:
        missing_display = ", ".join(sorted(fallback_failures))
        raise Auto3DExecutionError(
            "Auto3D returned optimized seeds for "
            f"{len(observed_ids)}/{len(protomer_records)} protomers; "
            f"missing outputs for: {missing_display}. "
            f"RDKit fallback also failed. Inspect Auto3D artifacts under {output_dir} "
            "(notably Auto3D.log and job*/smi_taut_3d.sdf)."
        )

    return records + fallback_records


def _rdkit_fallback_seed_for_missing_protomer(
    protomer: ProtomerRecord,
    *,
    config: RunConfig,
    command: list[str],
    output_sdf: Path,
    output_dir: Path,
    reason: str,
) -> SeedConformerRecord | None:
    fallback_config = config.model_copy(
        update={
            "seeding": config.seeding.model_copy(
                update={
                    "method": "etkdg",
                    "rdkit_num_conformers": 1,
                }
            )
        }
    )
    stereo = StereoRecord(
        id=protomer.id,
        parent_id=protomer.parent_id,
        input_molecule_id=protomer.input_molecule_id,
        molname=protomer.molname,
        canonical_smiles=protomer.canonical_smiles,
        isomeric_smiles=protomer.isomeric_smiles,
        molecular_formula=protomer.molecular_formula,
        formal_charge=protomer.formal_charge,
        explicit_proton_count=protomer.explicit_proton_count,
        source_software=protomer.source_software,
        source_command=protomer.source_command,
        source_python_function=(
            "dsvr.chemistry.conformers_auto3d._rdkit_fallback_seed_for_missing_protomer"
        ),
        warnings=list(protomer.warnings),
        metadata=dict(protomer.metadata),
        stereo_index=1,
        rdkit_mol=Chem.Mol(protomer.rdkit_mol),
    )
    rdkit_records = generate_rdkit_seeds(stereo, fallback_config)
    successful = [
        record
        for record in rdkit_records
        if record.embedding_status != "failed" and record.rdkit_mol is not None
    ]
    if not successful:
        return None
    record = successful[0]
    warning = (
        "Auto3D produced no optimized representative for this protomer; "
        "falling back to a single RDKit ETKDG seed. "
        f"Reason: {reason}"
    )
    metadata = dict(record.metadata)
    metadata["auto3d_fallback"] = {
        "auto3d_command": " ".join(command),
        "auto3d_output_sdf": str(output_sdf),
        "auto3d_job_dir": str(output_sdf.parent),
        "missing_protomer_id": protomer.id,
        "reason": reason,
        "cheap_3d_energy_omitted": True,
    }
    output_paths = list(record.output_paths)
    for extra_path in (
        output_sdf,
        output_dir / "auto3d_protocol_seeds.sdf",
        output_dir / "auto3d_protocol_seeds.csv",
    ):
        if extra_path not in output_paths:
            output_paths.append(extra_path)
    warnings = list(record.warnings)
    warnings.append(warning)
    return record.model_copy(
        update={
            "warnings": warnings,
            "metadata": metadata,
            "output_paths": output_paths,
            "energy_kcal_mol": None,
        }
    )


def _auto3d_partial_failure_reason(output_sdf: Path, output_dir: Path) -> str:
    candidate_logs = [output_sdf.parent / "Auto3D.log"]
    log_root = output_dir / "logs"
    if log_root.exists():
        candidate_logs.extend(sorted(log_root.glob("*_auto3d/stdout.log"), reverse=True))
        candidate_logs.extend(sorted(log_root.glob("*_auto3d/combined.log"), reverse=True))
    for log_path in candidate_logs:
        if not log_path.exists():
            continue
        try:
            text = log_path.read_text(encoding="utf-8")
        except OSError:
            continue
        if "Dropped(Oscillating)" in text:
            return "Auto3D dropped all candidate structures as oscillating during optimization"
    return "Auto3D did not emit an optimized structure for this protomer"


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


def _auto3d_variant_key(
    record: SeedConformerRecord,
) -> tuple[str, str, str, int | None, int | None]:
    return (
        record.input_molecule_id,
        record.isomeric_smiles or record.canonical_smiles or record.id,
        record.molecular_formula or "",
        record.formal_charge,
        record.explicit_proton_count,
    )


def _plausibility_population_group_key(
    record: SeedConformerRecord,
    config: RunConfig,
) -> tuple[object, ...]:
    if config.thermo.population_scope == "same_formula":
        return (
            "same_formula",
            record.molecular_formula,
            record.explicit_proton_count,
        )
    if config.thermo.population_scope == "same_charge":
        return ("same_charge", record.formal_charge)
    return ("all_approximate", "all")


def _write_auto3d_representative_outputs(
    output_dir: Path,
    records: list[ThermoRecord],
) -> None:
    csv_path = output_dir / "auto3d_representative_scores.csv"
    jsonl_path = output_dir / "auto3d_representative_scores.jsonl"
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
        "plausibility_score_kcal_mol",
        "auto3d_energy_kcal_mol",
        "candidate_conformer_count",
        "score_reasons",
        "warnings",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            representative = record.metadata.get("auto3d_representative", {})
            svp_score = record.metadata.get("svp_score", {})
            writer.writerow(
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
                    "temperature_kelvin": record.temperature_kelvin,
                    "plausibility_score_kcal_mol": record.free_energy_kcal_mol,
                    "auto3d_energy_kcal_mol": representative.get("auto3d_energy_kcal_mol"),
                    "candidate_conformer_count": representative.get("candidate_conformer_count"),
                    "score_reasons": " | ".join(svp_score.get("reasons", [])),
                    "warnings": " | ".join(record.warnings),
                }
            )
    with jsonl_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(record.model_dump_json() + "\n")


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
