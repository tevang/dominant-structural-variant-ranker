from __future__ import annotations

import csv
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from dsvr.config import RunConfig
from dsvr.filtering.budgets import seed_budget
from dsvr.models import SeedConformerRecord, StereoRecord, make_seed_id


def rdkit_available() -> bool:
    try:
        import rdkit  # noqa: F401
    except ImportError:
        return False
    return True


def generate_rdkit_seeds(
    stereo_record: StereoRecord,
    config: RunConfig,
) -> list[SeedConformerRecord]:
    output_dir = config.output_dir / "seeding" / "rdkit"
    xyz_dir = output_dir / "xyz"
    output_dir.mkdir(parents=True, exist_ok=True)
    write_raw_xyz = config.disk.keep_raw_xyz or config.disk.cleanup_policy == "debug_all"
    if write_raw_xyz:
        xyz_dir.mkdir(parents=True, exist_ok=True)

    mol = Chem.AddHs(Chem.Mol(stereo_record.rdkit_mol))
    params = _etkdg_params(config)
    status = "success"
    try:
        conformer_ids = list(
            AllChem.EmbedMultipleConfs(
                mol,
                numConfs=config.seeding.rdkit_num_conformers,
                params=params,
            )
        )
    except Exception as exc:
        return [
            _failure_record(
                stereo_record,
                config,
                output_dir,
                f"RDKit ETKDG embedding raised {type(exc).__name__}: {exc}",
            )
        ]

    if not conformer_ids:
        return [
            _failure_record(
                stereo_record,
                config,
                output_dir,
                "RDKit ETKDG returned zero conformers",
            )
        ]

    candidates: list[tuple[int, Chem.Mol, str, bool | None, float | None]] = []
    for index, conformer_id in enumerate(conformer_ids, start=1):
        seed_mol = Chem.Mol(mol)
        _remove_other_conformers(seed_mol, conformer_id)
        seed_mol.GetConformer().SetId(0)
        ff_status, converged, energy = _minimize(seed_mol, config)
        candidates.append((index, seed_mol, ff_status, converged, energy))

    selected_indexes, selection_rows = _select_seed_candidates(candidates, config)
    records: list[SeedConformerRecord] = []
    for selected_rank, (index, seed_mol, ff_status, converged, energy) in enumerate(
        [candidate for candidate in candidates if candidate[0] in selected_indexes],
        start=1,
    ):
        canonical_smiles = Chem.MolToSmiles(seed_mol, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(seed_mol, canonical=True, isomericSmiles=True)
        output_paths = [
            output_dir / f"{stereo_record.id}_seeds.sdf",
            output_dir / f"{stereo_record.id}_seeds.csv",
            output_dir / f"{stereo_record.id}_seed_selection.csv",
        ]
        if write_raw_xyz:
            xyz_path = xyz_dir / f"{stereo_record.id}_c{index:02d}.xyz"
            _write_xyz(seed_mol, xyz_path, comment=f"{stereo_record.id} conformer {index}")
            output_paths.append(xyz_path)
        metadata = {
            "rdkit_embedding": {
                "method": "ETKDGv3" if hasattr(AllChem, "ETKDGv3") else "ETKDG",
                "random_seed": config.enumeration.stereo_random_seed,
                "requested_conformers": config.seeding.rdkit_num_conformers,
                "prune_rms_thresh": config.seeding.rdkit_prune_rms_thresh,
                "raw_xyz_written": write_raw_xyz,
                "selected_seed_rank": selected_rank,
                "selected_seed_count": len(selected_indexes),
            },
            "forcefield": {
                "configured": config.seeding.rdkit_forcefield,
                "status": ff_status,
                "converged": converged,
            },
        }
        records.append(
            SeedConformerRecord(
                id=make_seed_id(
                    stereo_record.id,
                    selected_rank,
                    canonical_smiles,
                    isomeric_smiles,
                    metadata,
                ),
                parent_id=stereo_record.id,
                input_molecule_id=stereo_record.input_molecule_id,
                molname=stereo_record.molname,
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                molecular_formula=_formula(seed_mol),
                formal_charge=Chem.GetFormalCharge(seed_mol),
                explicit_proton_count=_explicit_proton_count(seed_mol),
                source_software="rdkit",
                source_python_function="dsvr.chemistry.conformers_rdkit.generate_rdkit_seeds",
                output_paths=output_paths,
                warnings=[],
                metadata=metadata,
                conformer_index=selected_rank,
                energy_kcal_mol=energy,
                rdkit_mol=seed_mol,
                rdkit_conformer_id=0,
                forcefield=config.seeding.rdkit_forcefield,
                forcefield_status=ff_status,
                minimization_converged=converged,
                embedding_status=status,
            )
        )

    _write_seed_sdf(output_dir / f"{stereo_record.id}_seeds.sdf", records)
    _write_seed_csv(output_dir / f"{stereo_record.id}_seeds.csv", records)
    _write_seed_selection_csv(output_dir / f"{stereo_record.id}_seed_selection.csv", selection_rows)
    return records


def read_stereo_sdf(path: Path) -> list[StereoRecord]:
    supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
    records: list[StereoRecord] = []
    for index, molecule in enumerate(supplier, start=1):
        if molecule is None:
            continue
        stereo_id = _prop_or_default(molecule, "DSVR_STEREO_ID", f"stereo_{index:06d}")
        input_id = _prop_or_default(molecule, "DSVR_INPUT_ID", stereo_id)
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
        isomeric_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        records.append(
            StereoRecord(
                id=stereo_id,
                parent_id=_prop_or_default(molecule, "DSVR_PARENT_TAUTOMER_ID", input_id),
                input_molecule_id=input_id,
                molname=_prop_or_default(molecule, "DSVR_MOLNAME", molecule.GetProp("_Name")),
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                molecular_formula=_formula(molecule),
                formal_charge=Chem.GetFormalCharge(molecule),
                explicit_proton_count=_explicit_proton_count(molecule),
                source_software="sdf",
                source_python_function="dsvr.chemistry.conformers_rdkit.read_stereo_sdf",
                stereo_index=index,
                rdkit_mol=molecule,
            )
        )
    return records


def _etkdg_params(config: RunConfig) -> AllChem.EmbedParameters:
    params = AllChem.ETKDGv3() if hasattr(AllChem, "ETKDGv3") else AllChem.ETKDG()
    params.randomSeed = config.enumeration.stereo_random_seed
    params.pruneRmsThresh = config.seeding.rdkit_prune_rms_thresh
    return params


def _minimize(molecule: Chem.Mol, config: RunConfig) -> tuple[str, bool | None, float | None]:
    forcefield = config.seeding.rdkit_forcefield
    if forcefield == "none":
        return "disabled", None, None
    if forcefield == "mmff":
        props = AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant="MMFF94")
        if props is None:
            return "mmff_parameters_unavailable", None, None
        ff = AllChem.MMFFGetMoleculeForceField(molecule, props)
        if ff is None:
            return "mmff_forcefield_unavailable", None, None
        result = ff.Minimize()
        return "mmff_minimized", result == 0, float(ff.CalcEnergy())

    ff = AllChem.UFFGetMoleculeForceField(molecule)
    if ff is None:
        return "uff_forcefield_unavailable", None, None
    result = ff.Minimize()
    return "uff_minimized", result == 0, float(ff.CalcEnergy())


def _failure_record(
    stereo_record: StereoRecord,
    config: RunConfig,
    output_dir: Path,
    error: str,
) -> SeedConformerRecord:
    metadata = {
        "rdkit_embedding": {
            "method": "ETKDGv3" if hasattr(AllChem, "ETKDGv3") else "ETKDG",
            "random_seed": config.enumeration.stereo_random_seed,
            "requested_conformers": config.seeding.rdkit_num_conformers,
            "prune_rms_thresh": config.seeding.rdkit_prune_rms_thresh,
        },
        "failure": error,
    }
    canonical_smiles = stereo_record.canonical_smiles
    isomeric_smiles = stereo_record.isomeric_smiles
    record = SeedConformerRecord(
        id=make_seed_id(stereo_record.id, 0, canonical_smiles, isomeric_smiles, metadata),
        parent_id=stereo_record.id,
        input_molecule_id=stereo_record.input_molecule_id,
        molname=stereo_record.molname,
        canonical_smiles=canonical_smiles,
        isomeric_smiles=isomeric_smiles,
        molecular_formula=stereo_record.molecular_formula,
        formal_charge=stereo_record.formal_charge,
        explicit_proton_count=stereo_record.explicit_proton_count,
        source_software="rdkit",
        source_python_function="dsvr.chemistry.conformers_rdkit.generate_rdkit_seeds",
        output_paths=[output_dir / f"{stereo_record.id}_seeds.csv"],
        warnings=[error],
        metadata=metadata,
        conformer_index=0,
        forcefield=config.seeding.rdkit_forcefield,
        forcefield_status="not_run",
        minimization_converged=None,
        embedding_status="failed",
    )
    _write_seed_csv(output_dir / f"{stereo_record.id}_seeds.csv", [record])
    return record


def _remove_other_conformers(molecule: Chem.Mol, keep_conformer_id: int) -> None:
    for conformer in list(molecule.GetConformers()):
        if conformer.GetId() != keep_conformer_id:
            molecule.RemoveConformer(conformer.GetId())


def _select_seed_candidates(
    candidates: list[tuple[int, Chem.Mol, str, bool | None, float | None]],
    config: RunConfig,
) -> tuple[set[int], list[dict[str, object]]]:
    max_seeds = seed_budget(config) or config.variant_filtering.max_seeds_per_variant
    if max_seeds <= 0:
        max_seeds = 1
    ranked = sorted(
        candidates,
        key=lambda item: (float("inf") if item[4] is None else item[4], item[0]),
    )
    selected: list[tuple[int, Chem.Mol, str, bool | None, float | None]] = []
    rows: list[dict[str, object]] = []
    for candidate in ranked:
        index, molecule, ff_status, converged, energy = candidate
        if len(selected) >= max_seeds:
            rows.append(
                _seed_selection_row(
                    candidate,
                    False,
                    "rejected_over_max_seeds_per_variant",
                    None,
                )
            )
            continue
        min_rms = _minimum_rms(molecule, [item[1] for item in selected])
        diverse = min_rms is None or min_rms >= config.seeding.rdkit_prune_rms_thresh
        if diverse or not selected:
            selected.append(candidate)
            rows.append(
                _seed_selection_row(
                    candidate,
                    True,
                    "selected_low_energy_rms_diverse",
                    min_rms,
                )
            )
        else:
            rows.append(
                _seed_selection_row(
                    candidate,
                    False,
                    "rejected_rms_duplicate",
                    min_rms,
                )
            )
    if len(selected) < max_seeds:
        selected_ids = {item[0] for item in selected}
        for candidate in ranked:
            if len(selected) >= max_seeds:
                break
            if candidate[0] in selected_ids:
                continue
            selected.append(candidate)
            selected_ids.add(candidate[0])
            _replace_selection_row(
                rows,
                candidate[0],
                _seed_selection_row(
                    candidate,
                    True,
                    "selected_low_energy_fill_to_seed_budget",
                    _minimum_rms(
                        candidate[1],
                        [item[1] for item in selected if item[0] != candidate[0]],
                    ),
                ),
            )
    selected_indexes = {item[0] for item in selected}
    return selected_indexes, sorted(
        rows,
        key=lambda row: int(str(row["raw_conformer_index"])),
    )


def _seed_selection_row(
    candidate: tuple[int, Chem.Mol, str, bool | None, float | None],
    selected: bool,
    reason: str,
    min_rms_to_selected: float | None,
) -> dict[str, object]:
    index, _molecule, ff_status, converged, energy = candidate
    return {
        "raw_conformer_index": index,
        "selected": selected,
        "reason": reason,
        "energy_kcal_mol": energy,
        "forcefield_status": ff_status,
        "minimization_converged": converged,
        "min_rms_to_selected": min_rms_to_selected,
    }


def _replace_selection_row(
    rows: list[dict[str, object]],
    raw_conformer_index: int,
    replacement: dict[str, object],
) -> None:
    for index, row in enumerate(rows):
        if row["raw_conformer_index"] == raw_conformer_index:
            rows[index] = replacement
            return
    rows.append(replacement)


def _minimum_rms(molecule: Chem.Mol, selected: list[Chem.Mol]) -> float | None:
    if not selected:
        return None
    values = []
    for selected_mol in selected:
        try:
            values.append(float(AllChem.GetBestRMS(molecule, selected_mol)))
        except (RuntimeError, ValueError):
            continue
    return min(values) if values else None


def _write_seed_sdf(path: Path, records: list[SeedConformerRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        if record.rdkit_mol is None:
            continue
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
            "DSVR_EXPLICIT_PROTON_COUNT": str(record.explicit_proton_count),
            "DSVR_FORCEFIELD": record.forcefield or "",
            "DSVR_FORCEFIELD_STATUS": record.forcefield_status,
            "DSVR_MINIMIZATION_CONVERGED": str(record.minimization_converged),
            "DSVR_ENERGY_KCAL_MOL": (
                "" if record.energy_kcal_mol is None else str(record.energy_kcal_mol)
            ),
        }.items():
            mol.SetProp(key, value)
        writer.write(mol)
    writer.close()


def _write_seed_csv(path: Path, records: list[SeedConformerRecord]) -> None:
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
        "forcefield",
        "forcefield_status",
        "minimization_converged",
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


def _write_seed_selection_csv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = [
        "raw_conformer_index",
        "selected",
        "reason",
        "energy_kcal_mol",
        "forcefield_status",
        "minimization_converged",
        "min_rms_to_selected",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column) for column in columns})


def _write_xyz(molecule: Chem.Mol, path: Path, comment: str) -> None:
    conformer = molecule.GetConformer()
    lines = [str(molecule.GetNumAtoms()), comment]
    for atom in molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol()} {position.x:.10f} {position.y:.10f} {position.z:.10f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _formula(molecule: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(molecule)


def _explicit_proton_count(molecule: Chem.Mol) -> int:
    with_hs = Chem.AddHs(molecule)
    return sum(1 for atom in with_hs.GetAtoms() if atom.GetAtomicNum() == 1)


def _prop_or_default(molecule: Chem.Mol, key: str, default: str) -> str:
    return molecule.GetProp(key) if molecule.HasProp(key) else default
