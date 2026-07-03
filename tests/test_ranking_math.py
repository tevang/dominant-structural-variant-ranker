import math

import pytest

from dsvr.ranking.boltzmann import R_KCAL_MOL_K, boltzmann_weights


def test_boltzmann_population_synthetic_delta_g_values() -> None:
    temperature = 298.15
    delta_g = [0.0, 1.0, 2.0]

    populations = boltzmann_weights(delta_g, temperature)

    raw = [math.exp(-value / (R_KCAL_MOL_K * temperature)) for value in delta_g]
    total = sum(raw)
    expected = [value / total for value in raw]
    assert populations == pytest.approx(expected)
    assert sum(populations) == pytest.approx(1.0)
    assert populations[0] > populations[1] > populations[2]
