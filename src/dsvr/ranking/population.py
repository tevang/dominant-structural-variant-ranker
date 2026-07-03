from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Protocol

from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.models import (
    CrestConformerRecord,
    RankedVariantRecord,
    ThermoRecord,
    make_ranked_variant_id,
)
from dsvr.ranking.boltzmann import boltzmann_weights
from dsvr.workflow.provenance import write_jsonl

CROSS_PROTOMER_WARNING = (
    "Population comparison across different protonation states is approximate without "
    "micro-pKa/proton chemical-potential corrections."
)


RankableRecord = ThermoRecord | CrestConformerRecord


class ProtonationCorrectionProvider(Protocol):
    """Future extension point for pH-dependent microstate free-energy corrections."""

    def get_microstate_correction(
        self,
        record: RankableRecord,
        ph: float,
        solvent: str,
        temperature: float,
    ) -> float | None:
        """Return an additive correction in kcal/mol, or None if unavailable."""
        ...


def approximate_populations(relative_energies_kcal_mol: list[float]) -> list[float]:
    return boltzmann_weights(relative_energies_kcal_mol)


def compute_delta_g_and_populations(
    records: list[RankableRecord],
    config: RunConfig,
    correction_provider: ProtonationCorrectionProvider | None = None,
) -> list[RankedVariantRecord]:
    corrections = _microstate_corrections(records, config, correction_provider)
    scored = [record for record in records if _score(record, corrections) is not None]
    grouped: dict[tuple[object, ...], list[RankableRecord]] = defaultdict(list)
    for record in scored:
        grouped[_population_group_key(record, config)].append(record)

    ranked: list[RankedVariantRecord] = []
    for group_key, group_records in grouped.items():
        group_scores = [_score(record, corrections) for record in group_records]
        finite_scores = [score for score in group_scores if score is not None]
        if not finite_scores:
            continue
        minimum = min(finite_scores)
        relative = [
            float(_score(record, corrections) or 0.0) - minimum for record in group_records
        ]
        populations = boltzmann_weights(relative, config.chemistry.temperature_kelvin)
        mixed_formula_or_protons = _has_mixed_formula_or_protons(group_records)
        corrections_complete = _corrections_complete_for_mixed_group(
            group_records,
            corrections,
            mixed_formula_or_protons,
        )
        for rank, (record, delta_g, population) in enumerate(
            sorted(
                zip(group_records, relative, populations, strict=True),
                key=lambda item: (item[1], item[0].id),
            ),
            start=1,
        ):
            warnings = list(record.warnings)
            approximate = _population_is_approximate(
                config,
                mixed_formula_or_protons,
                corrections_complete,
            )
            if mixed_formula_or_protons and not corrections_complete:
                warnings.append(CROSS_PROTOMER_WARNING)
            metadata = {
                "ranking": {
                    "group_key": list(group_key),
                    "population_scope": config.thermo.population_scope,
                    "temperature_kelvin": config.chemistry.temperature_kelvin,
                    "population_assumption": _population_assumption(
                        config,
                        mixed_formula_or_protons,
                        corrections_complete,
                    ),
                    "mixed_formula_or_proton_count": mixed_formula_or_protons,
                    "microstate_correction_kcal_mol": corrections.get(record.id),
                    "corrected_score_kcal_mol": _score(record, corrections),
                    "source_record_id": record.id,
                    "source_stage": record.stage_name,
                    "source_workdir": record.metadata.get("crest", {}).get("workdir"),
                }
            }
            ranked.append(
                RankedVariantRecord(
                    id=make_ranked_variant_id(
                        record.id,
                        rank,
                        record.canonical_smiles,
                        record.isomeric_smiles,
                        metadata,
                    ),
                    parent_id=record.id,
                    input_molecule_id=record.input_molecule_id,
                    molname=record.molname,
                    canonical_smiles=record.canonical_smiles,
                    isomeric_smiles=record.isomeric_smiles,
                    molecular_formula=record.molecular_formula,
                    formal_charge=record.formal_charge,
                    explicit_proton_count=record.explicit_proton_count,
                    source_software="dsvr-ranking",
                    source_python_function="dsvr.ranking.population.compute_delta_g_and_populations",
                    warnings=sorted(set(warnings)),
                    metadata=metadata,
                    rank=rank,
                    score_kcal_mol=_score(record, corrections),
                    relative_free_energy_kcal_mol=delta_g,
                    boltzmann_population=population,
                    population_scope=config.thermo.population_scope,
                    approximate_population=approximate,
                )
            )

    return sorted(
        ranked,
        key=lambda item: (
            str(item.metadata["ranking"]["group_key"]),
            item.rank,
            item.id,
        ),
    )


def write_ranked_outputs(records: list[RankedVariantRecord], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_ranked_csv(output_dir / "ranked_variants.csv", records)
    _write_ranked_json(output_dir / "ranked_variants.json", records)
    _write_ranked_sdf(output_dir / "ranked_variants.sdf", records)
    _write_ranking_summary(output_dir / "ranking_summary.md", records)
    write_jsonl(output_dir / "ranking_provenance.jsonl", records)


def load_rankable_records(run_dir: Path) -> list[RankableRecord]:
    thermo_records = [
        ThermoRecord.model_validate(json.loads(path.read_text(encoding="utf-8")))
        for path in sorted(run_dir.glob("**/xtb_thermo.json"))
    ]
    if thermo_records:
        return thermo_records

    crest_records: list[CrestConformerRecord] = []
    for path in sorted(run_dir.glob("**/crest_provenance.jsonl")):
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                data = json.loads(line)
                if data.get("stage_name") == "crest_conformer":
                    crest_records.append(CrestConformerRecord.model_validate(data))
    return crest_records


def _microstate_corrections(
    records: list[RankableRecord],
    config: RunConfig,
    correction_provider: ProtonationCorrectionProvider | None,
) -> dict[str, float]:
    if correction_provider is None:
        return {}
    corrections: dict[str, float] = {}
    for record in records:
        correction = correction_provider.get_microstate_correction(
            record,
            config.chemistry.ph,
            config.chemistry.solvent,
            config.chemistry.temperature_kelvin,
        )
        if correction is not None:
            corrections[record.id] = correction
    return corrections


def _score(record: RankableRecord, corrections: dict[str, float] | None = None) -> float | None:
    if isinstance(record, ThermoRecord):
        score = record.free_energy_kcal_mol
    else:
        score = record.energy_kcal_mol
    if score is None:
        return None
    return score + (corrections or {}).get(record.id, 0.0)


def _population_group_key(record: RankableRecord, config: RunConfig) -> tuple[object, ...]:
    if config.thermo.population_scope == "same_formula":
        return (
            "same_formula",
            record.molecular_formula,
            record.explicit_proton_count,
        )
    if config.thermo.population_scope == "same_charge":
        return ("same_charge", record.formal_charge)
    return ("all_approximate", "all")


def _has_mixed_formula_or_protons(records: list[RankableRecord]) -> bool:
    keys = {(record.molecular_formula, record.explicit_proton_count) for record in records}
    return len(keys) > 1


def _corrections_complete_for_mixed_group(
    records: list[RankableRecord],
    corrections: dict[str, float],
    mixed_formula_or_protons: bool,
) -> bool:
    return mixed_formula_or_protons and all(record.id in corrections for record in records)


def _population_is_approximate(
    config: RunConfig,
    mixed_formula_or_protons: bool,
    corrections_complete: bool,
) -> bool:
    if config.thermo.population_scope == "all_approximate":
        return True
    return mixed_formula_or_protons and not corrections_complete


def _population_assumption(
    config: RunConfig,
    mixed_formula_or_protons: bool,
    corrections_complete: bool,
) -> str:
    if config.thermo.population_scope == "same_formula" and not mixed_formula_or_protons:
        return "Comparable within same formula and explicit proton-count group."
    if mixed_formula_or_protons and corrections_complete:
        return (
            "Compared across formula/proton-count differences using supplied microstate "
            "corrections."
        )
    return (
        "Approximate across protonation/protomer or formula/proton-count differences "
        "unless micro-pKa/proton chemical-potential corrections are available."
    )


def _write_ranked_csv(path: Path, records: list[RankedVariantRecord]) -> None:
    columns = [
        "rank",
        "id",
        "parent_id",
        "input_molecule_id",
        "molname",
        "canonical_smiles",
        "isomeric_smiles",
        "molecular_formula",
        "formal_charge",
        "explicit_proton_count",
        "score_kcal_mol",
        "relative_free_energy_kcal_mol",
        "boltzmann_population",
        "population_scope",
        "population_is_approximate",
        "approximate_population",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = record.model_dump(mode="json")
            row["population_is_approximate"] = record.approximate_population
            row["warnings"] = " | ".join(record.warnings)
            writer.writerow({column: row.get(column) for column in columns})


def _write_ranked_json(path: Path, records: list[RankedVariantRecord]) -> None:
    path.write_text(
        json.dumps([record.model_dump(mode="json") for record in records], indent=2) + "\n",
        encoding="utf-8",
    )


def _write_ranked_sdf(path: Path, records: list[RankedVariantRecord]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        mol = None
        if record.isomeric_smiles:
            mol = Chem.MolFromSmiles(record.isomeric_smiles)
        if mol is None and record.canonical_smiles:
            mol = Chem.MolFromSmiles(record.canonical_smiles)
        if mol is None:
            continue
        mol.SetProp("_Name", record.id)
        for key, value in {
            "DSVR_RANK": str(record.rank),
            "DSVR_PARENT_ID": record.parent_id or "",
            "DSVR_INPUT_ID": record.input_molecule_id,
            "DSVR_MOLNAME": record.molname,
            "DSVR_SCORE_KCAL_MOL": str(record.score_kcal_mol),
            "DSVR_RELATIVE_FREE_ENERGY_KCAL_MOL": str(record.relative_free_energy_kcal_mol),
            "DSVR_BOLTZMANN_POPULATION": str(record.boltzmann_population),
            "DSVR_POPULATION_SCOPE": record.population_scope,
            "DSVR_APPROXIMATE_POPULATION": str(record.approximate_population),
            "DSVR_WARNINGS": " | ".join(record.warnings),
        }.items():
            mol.SetProp(key, value)
        writer.write(mol)
    writer.close()


def _write_ranking_summary(path: Path, records: list[RankedVariantRecord]) -> None:
    lines = [
        "# DSVR Ranking Summary",
        "",
        "Population estimates are derived from relative free energies via Boltzmann weights.",
        CROSS_PROTOMER_WARNING,
        "",
        "| Rank | Molecule | ΔG kcal/mol | Population | Scope | Approximate |",
        "| --- | --- | ---: | ---: | --- | --- |",
    ]
    for record in records:
        lines.append(
            "| "
            f"{record.rank} | {record.molname} | "
            f"{_format_optional(record.relative_free_energy_kcal_mol)} | "
            f"{record.boltzmann_population if record.boltzmann_population is not None else ''} | "
            f"{record.population_scope} | {record.approximate_population} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _format_optional(value: float | None) -> str:
    return "" if value is None else str(value)
