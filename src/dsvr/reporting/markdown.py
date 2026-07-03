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
            "- Ranking uses relative CREST/xTB-derived free energies or downstream "
            "refined energies when available."
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


def _markdown_to_simple_html(markdown: str) -> str:
    body = "\n".join(f"<pre>{html.escape(line)}</pre>" for line in markdown.splitlines())
    return f"<!doctype html><html><body>{body}</body></html>\n"
