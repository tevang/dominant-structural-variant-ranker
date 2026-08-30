"""2D structural depictions rendered as SVG from SMILES.

SDF/3D structural files are intentionally excluded from rendering per the
capability spec; these helpers only produce 2D depictions from SMILES strings
using RDKit (a base dependency).
"""

from __future__ import annotations

import functools

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


@functools.lru_cache(maxsize=256)
def smiles_to_svg(smiles: str, width: int = 260, height: int = 200) -> str | None:
    """Return an SVG string for a molecule, or None if the SMILES is invalid."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
