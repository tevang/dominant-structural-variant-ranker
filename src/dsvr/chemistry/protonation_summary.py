"""Derived protonation-ensemble properties from Uni-Pka microstate free energies.

Uni-Pka predicts one pH-independent free energy (dG) per microspecies. This
module re-derives Boltzmann occupancies at any pH (the same reweighting the
Uni-Pka script uses) and computes the uncertainty-facing summary properties:
ensemble entropy, charge populations, isoelectric point, and nearest pKa
transition distance. Transitions are located on a configurable pH grid and
refined by linear interpolation, so precision is bounded by grid step.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

LN10 = math.log(10.0)
# Effective room-temperature offset the Uni-Pka model applies when mapping
# microstate free energies to pH-dependent weights (hard-coded upstream).
TRANSLATE_PH = 6.504894871171601
@dataclass(frozen=True)
class Microstate:
    smiles: str
    charge: int
    dg: float


@dataclass(frozen=True)
class ProtonationSummary:
    top_two_occupancy_gap: float | None
    occupancy_entropy: float | None
    charge_population: dict[int, float] = field(default_factory=dict)
    microstate_count: int = 0
    pka_nearest_distance: float | None = None
    pka_nearest_transition: float | None = None
    isoelectric_point: float | None = None
    warnings: list[str] = field(default_factory=list)


def compute_occupancies(microstates: list[Microstate], ph: float) -> list[float]:
    """Boltzmann occupancies of microstates at a given pH (numerically stable)."""

    if not microstates:
        return []
    weights = [
        math.exp(-_weight_exponent(state.dg, state.charge, ph)) for state in microstates
    ]
    total = sum(weights)
    if total <= 0:
        return [1.0 / len(weights)] * len(weights)
    return [weight / total for weight in weights]


def compute_protonation_summary(
    microstates: list[Microstate],
    *,
    working_ph: float,
    selected_forms: list[tuple[str, float]],
    ph_low: float = 0.0,
    ph_high: float = 14.0,
    ph_step: float = 0.05,
) -> ProtonationSummary:
    """Compute all Uni-Pka summary properties for one input molecule.

    ``selected_forms`` are the forms that survived selection (SMILES, occupancy
    at working pH); order is irrelevant.
    """

    warnings: list[str] = []
    if not microstates:
        return ProtonationSummary(
            top_two_occupancy_gap=None,
            occupancy_entropy=None,
            warnings=["no Uni-Pka microstates available; summary properties null"],
        )

    working_occupancies = compute_occupancies(microstates, working_ph)
    entropy = -sum(o * math.log(o) for o in working_occupancies if o > 0.0)
    charge_population = _charge_population(microstates, working_occupancies)

    if len(selected_forms) >= 2:
        # caller order is selection priority (keep_input_state first), not
        # occupancy order, so rank explicitly before taking the top-two gap
        top_two = sorted((occupancy for _, occupancy in selected_forms), reverse=True)[:2]
        gap: float | None = top_two[0] - top_two[1]
    elif selected_forms:
        gap = selected_forms[0][1]
    else:
        gap = None

    grid = _build_grid(ph_low, ph_high, working_ph, ph_step)
    net_charge_curve = [
        _net_charge(microstates, compute_occupancies(microstates, pH)) for pH in grid
    ]

    isoelectric = _first_crossing(grid, net_charge_curve, target=0.0)
    # Macroscopic pKa transitions sit where the population-weighted net charge
    # crosses a half-integer level k+0.5 (the 50/50 point between the charge-k
    # and charge-(k+1) families); the isoelectric point is the zero-charge
    # crossing and is reported separately.
    transitions: list[float] = [
        t
        for level in _half_integer_levels_between(net_charge_curve)
        if (t := _first_crossing(grid, net_charge_curve, target=level)) is not None
    ]
    if transitions:
        nearest = min(transitions, key=lambda t: abs(t - working_ph))
        nearest_distance: float | None = abs(nearest - working_ph)
    else:
        nearest, nearest_distance = None, None

    if isoelectric is None:
        warnings.append(
            "net charge never changes sign over the scanned pH range; isoelectric point null"
        )
    if nearest is None:
        warnings.append(
            "no macroscopic charge-level transition found in the scanned pH range; "
            "pKa distance null"
        )

    return ProtonationSummary(
        top_two_occupancy_gap=gap,
        occupancy_entropy=entropy,
        charge_population=charge_population,
        microstate_count=len(microstates),
        pka_nearest_distance=nearest_distance,
        pka_nearest_transition=nearest,
        isoelectric_point=isoelectric,
        warnings=warnings,
    )


def _weight_exponent(dg: float, charge: int, ph: float) -> float:
    """Uni-Pka ensemble weight exponent: dG - ln(10)*(TRANSLATE_PH - pH)*q."""

    return dg - LN10 * (TRANSLATE_PH - ph) * charge


def _charge_population(
    microstates: list[Microstate], occupancies: list[float]
) -> dict[int, float]:
    population: dict[int, float] = {}
    for state, occupancy in zip(microstates, occupancies, strict=True):
        population[state.charge] = population.get(state.charge, 0.0) + occupancy
    return population


def _net_charge(microstates: list[Microstate], occupancies: list[float]) -> float:
    return sum(state.charge * occupancy for state, occupancy in zip(microstates, occupancies, strict=True))


def _build_grid(low: float, high: float, extra: float, step: float) -> list[float]:
    values: list[float] = []
    pH = low
    while pH <= high + 1e-9:
        values.append(round(pH, 6))
        pH += step
    if extra not in values:
        values.append(extra)
    return sorted(values)


def _first_crossing(grid: list[float], curve: list[float], target: float) -> float | None:
    """First pH where the curve crosses ``target`` with a sign change.

    A curve that merely sits at the target (e.g. a permanently neutral
    molecule vs. target 0) does not count as a crossing, since no transition
    occurs there.
    """

    eps = 1e-9

    def side(value: float) -> int:
        delta = value - target
        if delta > eps:
            return 1
        if delta < -eps:
            return -1
        return 0

    previous_side = side(curve[0])
    previous_value = curve[0] - target
    for index in range(1, len(grid)):
        current_side = side(curve[index])
        if previous_side != 0 and current_side != 0 and previous_side != current_side:
            value = curve[index] - target
            frac = abs(previous_value) / (abs(previous_value) + abs(value))
            return grid[index - 1] + frac * (grid[index] - grid[index - 1])
        if current_side != 0:
            previous_side = current_side
        previous_value = curve[index] - target
    return None


def _half_integer_levels_between(curve: list[float]) -> list[float]:
    """Half-integer charge levels k + 0.5 spanned by the curve, where k+0.5 is
    the midpoint of the transition between charge families k and k+1."""

    eps = 1e-9
    low = math.floor(min(curve) - eps)
    high = math.ceil(max(curve) + eps)
    return [
        level + 0.5
        for level in range(low, high)
        if min(curve) + eps < level + 0.5 < max(curve) - eps
    ]
