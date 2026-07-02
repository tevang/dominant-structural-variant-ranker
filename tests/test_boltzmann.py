import pytest

from dsvr.ranking.boltzmann import boltzmann_weights


def test_boltzmann_weights_normalize() -> None:
    weights = boltzmann_weights([0.0, 1.0])

    assert sum(weights) == pytest.approx(1.0)
    assert weights[0] > weights[1]

