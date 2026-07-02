from __future__ import annotations

from dsvr.ranking.boltzmann import boltzmann_weights


def approximate_populations(relative_energies_kcal_mol: list[float]) -> list[float]:
    return boltzmann_weights(relative_energies_kcal_mol)

