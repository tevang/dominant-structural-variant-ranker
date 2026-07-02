from __future__ import annotations

from pathlib import Path


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

