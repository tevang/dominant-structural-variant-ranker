from __future__ import annotations

import math

R_KCAL_MOL_K = 0.00198720425864083


def boltzmann_weights(
    relative_energies_kcal_mol: list[float],
    temperature_k: float = 298.15,
) -> list[float]:
    if temperature_k <= 0:
        raise ValueError("temperature_k must be positive")
    weights = [
        math.exp(-energy / (R_KCAL_MOL_K * temperature_k))
        for energy in relative_energies_kcal_mol
    ]
    total = sum(weights)
    return [weight / total for weight in weights] if total else [0.0 for _ in weights]
