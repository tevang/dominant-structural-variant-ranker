from __future__ import annotations

from pathlib import Path

from dsvr.io.sdf import read_sdf
from dsvr.io.smiles import read_smiles
from dsvr.models import InputMolecule


def read_molecules(path: Path) -> list[InputMolecule]:
    suffix = path.suffix.lower()
    if suffix in {".smi", ".smiles", ".txt"}:
        return read_smiles(path)
    if suffix in {".sdf", ".sd"}:
        return read_sdf(path)
    raise ValueError(f"Unsupported input format: {path}")

