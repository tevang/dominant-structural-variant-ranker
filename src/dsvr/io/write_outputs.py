from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.models import RankedVariantRecord, VariantRecord

RANKED_VARIANT_COLUMNS = [
    "rank",
    "molname",
    "input_id",
    "protomer_id",
    "tautomer_id",
    "stereo_id",
    "best_conformer_id",
    "canonical_smiles",
    "isomeric_smiles",
    "formula",
    "formal_charge",
    "proton_count",
    "solvent",
    "solvent_model",
    "pH",
    "temperature_kelvin",
    "electronic_energy",
    "thermo_correction",
    "free_energy",
    "relative_free_energy_kcal_mol",
    "boltzmann_weight",
    "population",
    "population_scope",
    "population_is_approximate",
    "ranking_method",
    "source_software",
    "warnings",
    "paths_to_structures",
]

SDF_RANKED_PROPERTIES = [
    "DSVR_RANK",
    "DSVR_MOLNAME",
    "DSVR_PROTOMER_ID",
    "DSVR_TAUTOMER_ID",
    "DSVR_STEREO_ID",
    "DSVR_CONFORMER_ID",
    "DSVR_FORMAL_CHARGE",
    "DSVR_FORMULA",
    "DSVR_PH",
    "DSVR_SOLVENT",
    "DSVR_RELATIVE_G_KCAL_MOL",
    "DSVR_POPULATION",
    "DSVR_POPULATION_SCOPE",
    "DSVR_POPULATION_APPROXIMATE",
    "DSVR_RANKING_METHOD",
    "DSVR_WARNINGS",
]


def write_ranked_csv(path: Path, records: list[VariantRecord]) -> None:
    rows = [
        record.model_dump(mode="json") | {"rank": rank}
        for rank, record in enumerate(records, 1)
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def write_final_ranked_outputs(
    output_dir: Path,
    records: list[RankedVariantRecord],
    config: RunConfig,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = [ranked_variant_row(record, config) for record in records]
    pd.DataFrame(rows, columns=RANKED_VARIANT_COLUMNS).to_csv(
        output_dir / "ranked_variants.csv",
        index=False,
    )
    write_json(
        output_dir / "ranked_variants.json",
        [record.model_dump(mode="json") for record in records],
    )
    write_ranked_sdf(output_dir / "ranked_variants.sdf", records, config)


def ranked_variant_row(record: RankedVariantRecord, config: RunConfig) -> dict[str, Any]:
    lineage = _lineage(record)
    source_paths = _paths_to_structures(record)
    score = record.score_kcal_mol
    electronic_energy = record.metadata.get("xtb", {}).get("electronic_energy_kcal_mol")
    thermo_correction = None
    if electronic_energy is not None and score is not None:
        thermo_correction = score - electronic_energy
    return {
        "rank": record.rank,
        "molname": record.molname,
        "input_id": record.input_molecule_id,
        "protomer_id": lineage.get("protomer_id"),
        "tautomer_id": lineage.get("tautomer_id"),
        "stereo_id": lineage.get("stereo_id"),
        "best_conformer_id": lineage.get("best_conformer_id") or record.parent_id,
        "canonical_smiles": record.canonical_smiles,
        "isomeric_smiles": record.isomeric_smiles,
        "formula": record.molecular_formula,
        "formal_charge": record.formal_charge,
        "proton_count": record.explicit_proton_count,
        "solvent": config.chemistry.solvent,
        "solvent_model": config.chemistry.solvent_model,
        "pH": config.chemistry.ph,
        "temperature_kelvin": config.chemistry.temperature_kelvin,
        "electronic_energy": electronic_energy,
        "thermo_correction": thermo_correction,
        "free_energy": score,
        "relative_free_energy_kcal_mol": record.relative_free_energy_kcal_mol,
        "boltzmann_weight": record.boltzmann_population,
        "population": record.boltzmann_population,
        "population_scope": record.population_scope,
        "population_is_approximate": record.approximate_population,
        "ranking_method": record.source_python_function or record.source_software,
        "source_software": record.source_software,
        "warnings": " | ".join(record.warnings),
        "paths_to_structures": json.dumps(source_paths, sort_keys=True),
    }


def write_ranked_sdf(
    path: Path,
    records: list[RankedVariantRecord],
    config: RunConfig,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    for record in records:
        mol = None
        if record.isomeric_smiles:
            mol = Chem.MolFromSmiles(record.isomeric_smiles)
        if mol is None and record.canonical_smiles:
            mol = Chem.MolFromSmiles(record.canonical_smiles)
        if mol is None:
            continue
        props = ranked_sdf_properties(record, config)
        mol.SetProp("_Name", record.id)
        for key, value in props.items():
            mol.SetProp(key, "" if value is None else str(value))
        writer.write(mol)
    writer.close()


def ranked_sdf_properties(
    record: RankedVariantRecord,
    config: RunConfig,
) -> dict[str, Any]:
    lineage = _lineage(record)
    return {
        "DSVR_RANK": record.rank,
        "DSVR_MOLNAME": record.molname,
        "DSVR_PROTOMER_ID": lineage.get("protomer_id"),
        "DSVR_TAUTOMER_ID": lineage.get("tautomer_id"),
        "DSVR_STEREO_ID": lineage.get("stereo_id"),
        "DSVR_CONFORMER_ID": lineage.get("best_conformer_id") or record.parent_id,
        "DSVR_FORMAL_CHARGE": record.formal_charge,
        "DSVR_FORMULA": record.molecular_formula,
        "DSVR_PH": config.chemistry.ph,
        "DSVR_SOLVENT": config.chemistry.solvent,
        "DSVR_RELATIVE_G_KCAL_MOL": record.relative_free_energy_kcal_mol,
        "DSVR_POPULATION": record.boltzmann_population,
        "DSVR_POPULATION_SCOPE": record.population_scope,
        "DSVR_POPULATION_APPROXIMATE": record.approximate_population,
        "DSVR_RANKING_METHOD": record.source_python_function or record.source_software,
        "DSVR_WARNINGS": " | ".join(record.warnings),
    }


def _lineage(record: RankedVariantRecord) -> dict[str, str | None]:
    metadata = record.metadata
    ranking = metadata.get("ranking", {})
    censo = metadata.get("censo", {})
    qm = metadata.get("qm", {})
    return {
        "protomer_id": ranking.get("protomer_id") or _ancestor_id(record.parent_id, "p"),
        "tautomer_id": ranking.get("tautomer_id") or _ancestor_id(record.parent_id, "t"),
        "stereo_id": ranking.get("stereo_id") or _ancestor_id(record.parent_id, "s"),
        "best_conformer_id": (
            ranking.get("source_record_id")
            or censo.get("preserves_preliminary_ranking_id")
            or qm.get("preserves_preliminary_ranking_id")
            or record.parent_id
        ),
    }


def _paths_to_structures(record: RankedVariantRecord) -> list[str]:
    paths: list[str] = []
    for section in ("ranking", "censo", "qm"):
        workdir = record.metadata.get(section, {}).get("source_workdir") or record.metadata.get(
            section,
            {},
        ).get("workdir")
        if workdir:
            paths.append(str(workdir))
    paths.extend(str(path) for path in record.output_paths)
    return sorted(set(paths))


def _ancestor_id(value: str | None, marker: str) -> str | None:
    if value is None:
        return None
    token = f"_{marker}"
    if token not in value:
        return None
    parts = value.split("_")
    for index, part in enumerate(parts):
        if part.startswith(marker) and len(part) >= 2 and part[1:3].isdigit():
            return "_".join(parts[: index + 2]) if index + 1 < len(parts) else value
    return None
