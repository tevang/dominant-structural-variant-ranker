from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from dsvr.config import RunConfig
from dsvr.models import SeedConformerRecord, StereoRecord, make_seed_id
from dsvr.runners.auto3d_runner import Auto3DExecutionError, Auto3DUnavailableError, run_auto3d

FINAL_3D_WARNING = (
    "Final Auto3D energies are approximate gas-phase/neural-potential conformer energies; "
    "they are not solvated free energies unless an optional validation stage is run."
)


@dataclass(frozen=True)
class Final3DResult:
    records: list[SeedConformerRecord]
    used_fallback: bool
    warnings: list[str]


def generate_final_3d_variants(
    stereo_records: list[StereoRecord],
    config: RunConfig,
) -> Final3DResult:
    output_dir = config.output_dir / "final_3d"
    output_dir.mkdir(parents=True, exist_ok=True)
    if not stereo_records:
        _write_final_outputs(config.output_dir, [], config)
        return Final3DResult(records=[], used_fallback=False, warnings=[])

    input_sdf = output_dir / "final_3d_input.sdf"
    _write_final_auto3d_input(input_sdf, stereo_records)
    try:
        output_sdf, command = _run_final_auto3d(input_sdf, output_dir, config)
        records = _records_from_final_auto3d_output(output_sdf, stereo_records, config, command)
        records = _dedupe_one_conformer_per_variant(records)
        records, fallback_warnings = _fill_missing_with_rdkit(
            stereo_records,
            records,
            config=config,
            command=command,
            output_sdf=output_sdf,
            reason="Auto3D returned no final conformer for this selected variant",
        )
        used_fallback = bool(fallback_warnings)
        warnings = fallback_warnings
    except (Auto3DExecutionError, Auto3DUnavailableError) as exc:
        command = ["auto3d", "unavailable_or_failed", str(input_sdf)]
        records = [_rdkit_final_record(record, config, command, input_sdf, str(exc)) for record in stereo_records]
        used_fallback = True
        warnings = [f"Auto3D final 3D generation failed; used RDKit one-conformer fallback: {exc}"]

    _write_final_outputs(config.output_dir, records, config)
    return Final3DResult(records=records, used_fallback=used_fallback, warnings=warnings)


def _run_final_auto3d(
    input_sdf: Path,
    output_dir: Path,
    config: RunConfig,
) -> tuple[Path, list[str]]:
    try:
        return run_auto3d(
            input_sdf,
            output_dir,
            k=config.final_3d.k,
            model=config.final_3d.optimizing_engine,
            internal_tautomer_stereo_enum=False,
            max_confs=config.final_3d.max_confs,
            patience=config.final_3d.patience,
            use_gpu=config.final_3d.use_gpu,
            stream_output=False,
            timeout_s=config.final_3d.timeout_seconds_per_batch,
        )
    except (Auto3DExecutionError, Auto3DUnavailableError):
        fallback = config.final_3d.fallback_optimizing_engine
        if fallback == config.final_3d.optimizing_engine:
            raise
        return run_auto3d(
            input_sdf,
            output_dir,
            k=config.final_3d.k,
            model=fallback,
            internal_tautomer_stereo_enum=False,
            max_confs=config.final_3d.max_confs,
            patience=config.final_3d.patience,
            use_gpu=config.final_3d.use_gpu,
            stream_output=False,
            timeout_s=config.final_3d.timeout_seconds_per_batch,
        )


def _write_final_auto3d_input(path: Path, records: list[StereoRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
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


def _records_from_final_auto3d_output(
    output_sdf: Path,
    stereo_records: list[StereoRecord],
    config: RunConfig,
    command: list[str],
) -> list[SeedConformerRecord]:
    stereo_by_id = {record.id: record for record in stereo_records}
    fallback = stereo_records[0] if len(stereo_records) == 1 else None
    supplier = Chem.SDMolSupplier(str(output_sdf), sanitize=True, removeHs=False)
    records: list[SeedConformerRecord] = []
    for index, molecule in enumerate(supplier, start=1):
        if molecule is None:
            continue
        stereo_id = _source_stereo_id(molecule)
        parent = stereo_by_id.get(stereo_id or "") or fallback
        if parent is None:
            continue
        energy, energy_prop = _extract_energy(molecule)
        records.append(
            _final_record_from_mol(
                molecule,
                parent,
                config,
                command,
                output_sdf,
                index=index,
                energy=energy,
                energy_property=energy_prop,
                fallback_reason=None,
            )
        )
    return records


def _final_record_from_mol(
    molecule: Chem.Mol,
    parent: StereoRecord,
    config: RunConfig,
    command: list[str],
    source_sdf: Path,
    *,
    index: int,
    energy: float | None,
    energy_property: str | None,
    fallback_reason: str | None,
) -> SeedConformerRecord:
    canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    isomeric_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    metadata = {
        "final_3d": {
            "tool": config.final_3d.tool,
            "k": config.final_3d.k,
            "max_confs": config.final_3d.max_confs,
            "patience": config.final_3d.patience,
            "optimizing_engine": config.final_3d.optimizing_engine,
            "fallback_optimizing_engine": config.final_3d.fallback_optimizing_engine,
            "use_gpu": config.final_3d.use_gpu,
            "one_conformer_per_variant": config.final_3d.one_conformer_per_variant,
            "energy_property": energy_property,
            "energy_warning": FINAL_3D_WARNING,
            "fallback_reason": fallback_reason,
        }
    }
    warnings = [FINAL_3D_WARNING]
    if fallback_reason:
        warnings.append(f"RDKit fallback final 3D conformer used: {fallback_reason}")
    return SeedConformerRecord(
        id=make_seed_id(parent.id, 1, canonical_smiles, isomeric_smiles, metadata),
        parent_id=parent.id,
        input_molecule_id=parent.input_molecule_id,
        molname=parent.molname,
        canonical_smiles=canonical_smiles,
        isomeric_smiles=isomeric_smiles,
        molecular_formula=_formula(molecule),
        formal_charge=Chem.GetFormalCharge(molecule),
        explicit_proton_count=_explicit_proton_count(molecule),
        source_software="rdkit" if fallback_reason else "auto3d",
        source_command=" ".join(command),
        source_python_function="dsvr.chemistry.final3d.generate_final_3d_variants",
        output_paths=[source_sdf, config.output_dir / "final_variants.sdf"],
        warnings=warnings,
        metadata=metadata,
        conformer_index=1,
        energy_kcal_mol=energy,
        rdkit_mol=molecule,
        rdkit_conformer_id=0 if molecule.GetNumConformers() else None,
        forcefield="rdkit" if fallback_reason else "auto3d",
        forcefield_status="rdkit_fallback" if fallback_reason else "auto3d_optimized",
        minimization_converged=None,
        embedding_status="success",
    )


def _dedupe_one_conformer_per_variant(records: list[SeedConformerRecord]) -> list[SeedConformerRecord]:
    best: dict[str, SeedConformerRecord] = {}
    for record in records:
        parent_id = record.parent_id or record.id
        current = best.get(parent_id)
        if current is None or _energy_sort_key(record) < _energy_sort_key(current):
            best[parent_id] = record
    return [best[key] for key in sorted(best)]


def _fill_missing_with_rdkit(
    stereo_records: list[StereoRecord],
    records: list[SeedConformerRecord],
    *,
    config: RunConfig,
    command: list[str],
    output_sdf: Path,
    reason: str,
) -> tuple[list[SeedConformerRecord], list[str]]:
    observed = {record.parent_id for record in records if record.parent_id}
    missing = [record for record in stereo_records if record.id not in observed]
    fallback_records = [
        _rdkit_final_record(record, config, command, output_sdf, reason) for record in missing
    ]
    warnings = [
        f"Auto3D produced no final conformer for {record.id}; used RDKit fallback."
        for record in missing
    ]
    return records + fallback_records, warnings


def _rdkit_final_record(
    parent: StereoRecord,
    config: RunConfig,
    command: list[str],
    source_sdf: Path,
    reason: str,
) -> SeedConformerRecord:
    mol = Chem.AddHs(Chem.Mol(parent.rdkit_mol))
    params = AllChem.ETKDGv3() if hasattr(AllChem, "ETKDGv3") else AllChem.ETKDG()
    params.randomSeed = int(config.enumeration.stereo_random_seed)
    conformer_id = int(AllChem.EmbedMolecule(mol, params))
    if conformer_id >= 0:
        mol.GetConformer(conformer_id).SetId(0)
    else:
        AllChem.Compute2DCoords(mol)
    energy = _minimize_rdkit(mol)
    return _final_record_from_mol(
        mol,
        parent,
        config,
        command,
        source_sdf,
        index=1,
        energy=energy,
        energy_property="rdkit_fallback_forcefield",
        fallback_reason=reason,
    )


def _write_final_outputs(output_dir: Path, records: list[SeedConformerRecord], config: RunConfig) -> None:
    _write_final_sdf(output_dir / "final_variants.sdf", records, config)
    _write_final_csv(output_dir / "final_variants.csv", records, config)
    _write_final_json(output_dir / "final_variants.json", records, config)
    _write_final_energy_csv(output_dir / "final_variant_energies.csv", records)


def _write_final_sdf(path: Path, records: list[SeedConformerRecord], config: RunConfig) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        mol = Chem.Mol(record.rdkit_mol)
        mol.SetProp("_Name", record.id)
        for key, value in _final_properties(record, config).items():
            mol.SetProp(key, "" if value is None else str(value))
        writer.write(mol)
    writer.close()


def _write_final_csv(path: Path, records: list[SeedConformerRecord], config: RunConfig) -> None:
    columns = [
        "final_variant_id",
        "molname",
        "input_id",
        "protomer_id",
        "tautomer_id",
        "stereoisomer_id",
        "canonical_smiles",
        "isomeric_smiles",
        "final_auto3d_energy_kcal_mol",
        "approximate_ranking",
        "warnings",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            props = _final_properties(record, config)
            writer.writerow(
                {
                    "final_variant_id": record.id,
                    "molname": record.molname,
                    "input_id": record.input_molecule_id,
                    "protomer_id": props["DSVR_PROTOMER_ID"],
                    "tautomer_id": props["DSVR_TAUTOMER_ID"],
                    "stereoisomer_id": props["DSVR_STEREO_ID"],
                    "canonical_smiles": record.canonical_smiles,
                    "isomeric_smiles": record.isomeric_smiles,
                    "final_auto3d_energy_kcal_mol": record.energy_kcal_mol,
                    "approximate_ranking": True,
                    "warnings": " | ".join(record.warnings),
                }
            )


def _write_final_json(path: Path, records: list[SeedConformerRecord], config: RunConfig) -> None:
    payload = [
        {
            "final_variant_id": record.id,
            "properties": _final_properties(record, config),
            "record": record.model_dump(mode="json"),
        }
        for record in records
    ]
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def _write_final_energy_csv(path: Path, records: list[SeedConformerRecord]) -> None:
    columns = ["final_variant_id", "stereoisomer_id", "energy_kcal_mol", "source_software"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "final_variant_id": record.id,
                    "stereoisomer_id": record.parent_id,
                    "energy_kcal_mol": record.energy_kcal_mol,
                    "source_software": record.source_software,
                }
            )


def _final_properties(record: SeedConformerRecord, config: RunConfig) -> dict[str, object | None]:
    protomer_id = _ancestor_id(record.parent_id, "p") or record.parent_id
    tautomer_id = _ancestor_id(record.parent_id, "t") or record.parent_id
    stereo_id = _ancestor_id(record.parent_id, "s") or record.parent_id
    return {
        "DSVR_INPUT_MOLNAME": record.molname,
        "DSVR_INPUT_ID": record.input_molecule_id,
        "DSVR_PROTOMER_ID": protomer_id,
        "DSVR_TAUTOMER_ID": tautomer_id,
        "DSVR_STEREO_ID": stereo_id,
        "DSVR_FINAL_VARIANT_ID": record.id,
        "DSVR_PROTONATION_STAGE": "molscrub_or_input_state",
        "DSVR_TAUTOMER_RELATIVE_ENERGY_KCAL_MOL": _metadata_energy(record, "tautomer_relative_energy_kcal_mol"),
        "DSVR_STEREO_RELATIVE_ENERGY_KCAL_MOL": _metadata_energy(record, "stereo_relative_energy_kcal_mol"),
        "DSVR_FINAL_AUTO3D_ENERGY_KCAL_MOL": record.energy_kcal_mol,
        "DSVR_APPROXIMATE_RANKING": True,
        "DSVR_ENERGY_WARNING": FINAL_3D_WARNING,
        "DSVR_WARNINGS": " | ".join(record.warnings),
        "DSVR_PH": config.chemistry.ph,
        "DSVR_SOLVENT": config.chemistry.solvent,
        "DSVR_SOLVENT_MODEL": config.chemistry.solvent_model,
    }


def _extract_energy(molecule: Chem.Mol) -> tuple[float | None, str | None]:
    for prop in (
        "E_kcal_mol",
        "E_rel(kcal/mol)",
        "E_tot(kcal/mol)",
        "energy_kcal_mol",
        "Energy",
        "ENERGY",
    ):
        if molecule.HasProp(prop):
            try:
                return float(molecule.GetProp(prop)), prop
            except ValueError:
                continue
    return None, None


def _source_stereo_id(molecule: Chem.Mol) -> str | None:
    for prop in ("DSVR_STEREO_ID", "DSVR_PARENT_STEREO_ID", "stereo_id", "parent_id", "_Name"):
        if molecule.HasProp(prop):
            value = molecule.GetProp(prop).strip()
            if value:
                return value
    return None


def _metadata_energy(record: SeedConformerRecord, key: str) -> float | None:
    for section in ("tautomer_filtering", "stereo_energy_filtering", "final_3d"):
        value = record.metadata.get(section, {}).get(key)
        if isinstance(value, int | float):
            return float(value)
    return None


def _energy_sort_key(record: SeedConformerRecord) -> tuple[float, str]:
    return (float("inf") if record.energy_kcal_mol is None else record.energy_kcal_mol, record.id)


def _minimize_rdkit(molecule: Chem.Mol) -> float | None:
    try:
        ff = AllChem.UFFGetMoleculeForceField(molecule)
        if ff is None:
            return None
        ff.Minimize()
        return float(ff.CalcEnergy())
    except (RuntimeError, ValueError):
        return None


def _formula(molecule: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(molecule)


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    return sum(atom.GetTotalNumHs(includeNeighbors=True) for atom in molecule.GetAtoms())


def _ancestor_id(value: str | None, marker: str) -> str | None:
    if value is None:
        return None
    parts = value.split("_")
    for index, part in enumerate(parts):
        if part.startswith(marker) and len(part) >= 2 and part[1:3].isdigit():
            return "_".join(parts[: index + 1])
    return None
