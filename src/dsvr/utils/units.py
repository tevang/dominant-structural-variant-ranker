from __future__ import annotations

HARTREE_TO_KCAL_MOL = 627.509474


def hartree_to_kcal_mol(value: float) -> float:
    return value * HARTREE_TO_KCAL_MOL

