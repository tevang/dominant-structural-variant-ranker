from __future__ import annotations

from pathlib import Path


def has_completed_smoke_outputs(outdir: Path) -> bool:
    return (outdir / "ranked.csv").exists() and (outdir / "provenance.json").exists()

