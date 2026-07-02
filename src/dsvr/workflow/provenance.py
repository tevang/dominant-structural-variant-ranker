from __future__ import annotations

import platform
import sys
from pathlib import Path

from dsvr.config import DsvrConfig
from dsvr.models import InputMolecule


def build_provenance(
    input_path: Path,
    config: DsvrConfig,
    molecules: list[InputMolecule],
) -> dict[str, object]:
    return {
        "input_path": str(input_path),
        "config": config.model_dump(mode="json"),
        "molecule_count": len(molecules),
        "input_hashes": [mol.input_hash for mol in molecules],
        "python": sys.version,
        "platform": platform.platform(),
        "scientific_limitation": (
            "pH controls candidate generation by default; cross-protomer "
            "populations are approximate "
            "without explicit micro-pKa/proton chemical-potential corrections."
        ),
    }
