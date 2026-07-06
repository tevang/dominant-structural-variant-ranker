from __future__ import annotations

import html
from collections import Counter
from pathlib import Path
from typing import Any

from dsvr.config import RunConfig
from dsvr.io.write_outputs import ranked_variant_row
from dsvr.models import AnyLineageRecord, RankedVariantRecord
from dsvr.utils.tool_check import check_tools


def write_summary_markdown(path: Path, molecule_count: int, variant_count: int) -> None:
    path.write_text(
        "\n".join(
            [
                "# DSVR Summary",
                "",
                f"- Molecules: {molecule_count}",
                f"- Variants: {variant_count}",
                "- Population estimates: approximate over generated candidates only.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_run_report(
    path: Path,
    *,
    config: RunConfig,
    records: list[AnyLineageRecord],
    ranked_records: list[RankedVariantRecord],
    manifest: dict[str, Any],
    output_files: list[Path],
    html_path: Path | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = _report_text(
        config=config,
        records=records,
        ranked_records=ranked_records,
        manifest=manifest,
        output_files=output_files,
    )
    path.write_text(text, encoding="utf-8")
    if html_path is not None:
        html_path.write_text(_markdown_to_simple_html(text), encoding="utf-8")


def _report_text(
    *,
    config: RunConfig,
    records: list[AnyLineageRecord],
    ranked_records: list[RankedVariantRecord],
    manifest: dict[str, Any],
    output_files: list[Path],
) -> str:
    counts = Counter(record.stage_name for record in records)
    tool_rows = [
        f"- {status.name} ({status.kind}): {'ok' if status.available else 'missing'} "
        f"{status.version or ''} {status.detail}"
        for status in check_tools(output_dir=config.output_dir)
    ]
    top_rows = []
    for record in sorted(ranked_records, key=lambda item: (item.input_molecule_id, item.rank)):
        if record.rank > 3:
            continue
        row = ranked_variant_row(record, config)
        top_rows.append(
            "| "
            f"{row['input_id']} | {row['rank']} | {row['molname']} | "
            f"{row['relative_free_energy_kcal_mol']} | {row['population']} | "
            f"{row['population_scope']} |"
        )
    warnings = sorted(
        {
            warning
            for record in records
            for warning in getattr(record, "warnings", [])
        }
        | set(manifest.get("warnings", []))
    )
    failures = [
        warning
        for warning in warnings
        if any(token in warning.lower() for token in ("failed", "unavailable", "error"))
    ]
    tautomer_timeout_rows = _tautomer_timeout_rows(records)
    stereo_energy = manifest.get("filtering", {}).get("stereo_energy_filtering", {})
    stereo_reduction = manifest.get("filtering", {}).get("stereo_reduction", {})
    xtb_prefilter = manifest.get("filtering", {}).get("xtb_prefilter", {})
    lines = [
        "# DSVR Run Report",
        "",
        "## Run Settings",
        "",
        f"- Run name: {config.run_name}",
        f"- Input: {config.input_path}",
        f"- Output directory: {config.output_dir}",
        f"- pH: {config.chemistry.ph}",
        f"- Solvent: {config.chemistry.solvent}",
        f"- Solvent model: {config.chemistry.solvent_model}",
        f"- Temperature K: {config.chemistry.temperature_kelvin}",
        f"- Population scope: {config.thermo.population_scope}",
        f"- Thermo max variants per molecule: {config.thermo.max_variants_per_molecule}",
        f"- Thermo max conformers per variant: {config.thermo.max_conformers_per_variant}",
        f"- CENSO enabled: {config.refinement.censo_enabled}",
        f"- CENSO max candidates: {config.refinement.max_candidates_for_refinement}",
        f"- QM enabled: {config.qm.enabled}",
        f"- QM max candidates: {config.qm.max_candidates}",
        f"- Final 3D tool: {config.final_3d.tool}",
        f"- Final 3D k: {config.final_3d.k}",
        f"- Final 3D max conformers: {config.final_3d.max_confs}",
        f"- Final 3D timeout seconds per batch: {config.final_3d.timeout_seconds_per_batch}",
        (
            "- RDKit tautomer timeout seconds: "
            f"{config.tautomer_filtering.rdkit_tautomer_timeout_seconds}"
        ),
        (
            "- RDKit tautomer max candidates before Auto3D: "
            f"{config.tautomer_filtering.max_rdkit_tautomers_before_auto3d}"
        ),
        f"- RDKit tautomer max transforms: {config.tautomer_filtering.max_rdkit_transforms}",
        f"- RDKit tautomer timeout fallback: {config.tautomer_filtering.fallback_if_timeout}",
        (
            "- RDKit stereoisomer timeout seconds: "
            f"{config.stereoisomer_filtering.timeout_seconds_per_tautomer}"
        ),
        (
            "- RDKit stereoisomer max isomers per tautomer: "
            f"{config.stereoisomer_filtering.max_stereoisomers_per_tautomer}"
        ),
        (
            "- Stereo Auto3D energy window kcal/mol: "
            f"{config.stereoisomer_filtering.stereo_energy_window_kcal_mol}"
        ),
        (
            "- Stereo keep top N diastereomers: "
            f"{config.stereoisomer_filtering.keep_top_n_diastereomers}"
        ),
        "",
        "## Tool Versions",
        "",
        *tool_rows,
        "",
        "## Stage Counts",
        "",
        *[f"- {stage}: {count}" for stage, count in sorted(counts.items())],
        "",
        "## Top Ranked Variants Per Input",
        "",
        "| Input ID | Rank | Molecule | ΔG kcal/mol | Population | Scope |",
        "| --- | ---: | --- | ---: | ---: | --- |",
        *(top_rows or ["| none | | | | | |"]),
        "",
        "## Scientific Assumptions and Warnings",
        "",
        (
            f"- pH {config.chemistry.ph} is used for candidate generation by default; "
            "molscrub-generated protonation/protomer states are not assigned rigorous "
            "pH populations."
        ),
        (
            f"- Solvent '{config.chemistry.solvent}' with solvent model "
            f"'{config.chemistry.solvent_model}' is used for the configured energy model."
        ),
        (
            "- Ranking uses approximate final Auto3D conformer energies in the default "
            "LigPrep-like path; optional CREST/xTB or downstream refined energies are "
            "used only when validation stages are enabled."
        ),
        (
            f"- Population scope is '{config.thermo.population_scope}'; states outside "
            "that scope are grouped separately or explicitly marked approximate."
        ),
        (
            "- Population comparison across different protonation states is approximate "
            "without micro-pKa/proton chemical-potential corrections."
        ),
        *[f"- {warning}" for warning in warnings],
        "",
        "## Tautomer Timeout Counts",
        "",
        "| Input ID | Molecule | Timeout fallback count |",
        "| --- | --- | ---: |",
        *(tautomer_timeout_rows or ["| none | none | 0 |"]),
        "",
        "## Stereoisomer Energy Filtering",
        "",
        f"- Enumerated stereo states: {stereo_energy.get('enumerated_count', 0)}",
        f"- Selected stereo states: {stereo_energy.get('selected_count', 0)}",
        f"- Rejected stereo states: {stereo_energy.get('rejected_count', 0)}",
        (
            "- Enantiomer states collapsed for Auto3D energy evaluation: "
            f"{stereo_energy.get('collapsed_count', 0)}"
        ),
        (
            "- Auto3D stereo energy evaluations run: "
            f"{stereo_energy.get('energy_evaluation_count', 0)}"
        ),
        "",
        "## Stereo Reduction",
        "",
        (
            "- CREST jobs saved by enantiomer collapse: "
            f"{stereo_reduction.get('jobs_saved', 0)}"
        ),
        (
            "- Enantiomer collapse assumes an achiral solvent/environment; disable "
            "`stereo_filtering.collapse_enantiomers_in_achiral_solvent` or set "
            "`stereo_filtering.solvent_is_chiral` for chiral binding-pocket evaluations."
        ),
        "",
        "## xTB Prefilter",
        "",
        f"- Enabled: {xtb_prefilter.get('enabled', False)}",
        f"- Decisions: {xtb_prefilter.get('decision_count', 0)}",
        f"- Variants/seeds pruned before CREST: {xtb_prefilter.get('pruned_count', 0)}",
        "",
        "## Failure Summary",
        "",
        *([f"- {failure}" for failure in failures] or ["- No failures recorded."]),
        "",
        "## Output File Locations",
        "",
        *[f"- {path}" for path in output_files],
        "",
    ]
    return "\n".join(lines)


def _tautomer_timeout_rows(records: list[AnyLineageRecord]) -> list[str]:
    counts: Counter[tuple[str, str]] = Counter()
    for record in records:
        if record.stage_name != "tautomer":
            continue
        if any(
            "tautomer enumeration timeout" in warning.lower()
            or "TAUTOMER_TIMEOUT_FALLBACK" in warning
            for warning in record.warnings
        ):
            counts[(record.input_molecule_id, record.molname)] += 1
    return [
        f"| {input_id} | {molname} | {count} |"
        for (input_id, molname), count in sorted(counts.items())
    ]


def _markdown_to_simple_html(markdown: str) -> str:
    body = "\n".join(f"<pre>{html.escape(line)}</pre>" for line in markdown.splitlines())
    return f"<!doctype html><html><body>{body}</body></html>\n"
