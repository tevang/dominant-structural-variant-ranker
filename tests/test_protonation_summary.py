"""Unit tests for Uni-Pka protonation summary math.

The reference formula (verified against the EasyDock Uni-Pka implementation's
documented imidazole example) is:

    weight(state, pH) ∝ exp(-dG - charge*LN10*(pH - TRANSLATE_PH))

so a pair (charged dG_c, neutral dG=0) is 50/50 at
pH = TRANSLATE_PH - dG_c/(charge*LN10).
"""

from __future__ import annotations

import math

import pytest

from dsvr.chemistry.protonation_summary import (
    LN10,
    TRANSLATE_PH,
    Microstate,
    compute_occupancies,
    compute_protonation_summary,
)


def dg_for_transition(ph_target: float, q: int) -> float:
    """dG of species with charge ``q`` at 50/50 vs neutral (dG=0) at ``ph_target``."""
    return -q * LN10 * (ph_target - TRANSLATE_PH)


def test_imidazole_documented_example() -> None:
    # EasyDock docs: dG=-6.7938 (q=+1), dG=-5.24537 (q=0); at pH 7.4 occupancies
    # are 0.3746 / 0.6254.
    states = [
        Microstate(smiles="prot", charge=1, dg=-6.7938),
        Microstate(smiles="neut", charge=0, dg=-5.24537),
    ]
    occ = compute_occupancies(states, ph=7.4)
    assert occ[0] == pytest.approx(0.3746, abs=1e-3)
    assert occ[1] == pytest.approx(0.6254, abs=1e-3)


def test_single_form_has_entropy_zero_and_gap_null() -> None:
    states = [Microstate(smiles="CCO", charge=0, dg=-3.0)]
    summary = compute_protonation_summary(states, working_ph=7.0, selected_forms=[("CCO", 1.0)])
    assert summary.occupancy_entropy == pytest.approx(0.0)
    assert summary.top_two_occupancy_gap == pytest.approx(1.0)
    assert summary.charge_population == {0: 1.0}
    assert summary.microstate_count == 1
    assert any("isoelectric" in w for w in summary.warnings)


def test_fifty_fifty_tie_entropy_is_ln2_and_gap_zero() -> None:
    states = [
        Microstate(smiles="acc", charge=0, dg=0.0),
        Microstate(smiles="accb", charge=0, dg=0.0),
    ]
    summary = compute_protonation_summary(
        states, working_ph=7.0, selected_forms=[("acc", 0.5), ("accb", 0.5)]
    )
    assert summary.occupancy_entropy == pytest.approx(math.log(2.0))
    assert summary.top_two_occupancy_gap == pytest.approx(0.0)


def test_acid_transition_located_near_pka() -> None:
    # Acid pair: neutral dG=0, anion with 50/50 at pH 4.76.
    states = [
        Microstate(smiles="A", charge=-1, dg=dg_for_transition(4.76, -1)),
        Microstate(smiles="B", charge=0, dg=0.0),
    ]
    occ = compute_occupancies(states, ph=4.76)
    assert occ[0] == pytest.approx(0.5, abs=1e-6)
    summary = compute_protonation_summary(
        states, working_ph=7.0, selected_forms=[("A", 0.999), ("B", 0.001)], ph_step=0.01
    )
    assert summary.pka_nearest_transition == pytest.approx(4.76, abs=0.01)
    assert summary.pka_nearest_distance == pytest.approx(7.0 - 4.76, abs=0.01)
    # A pure acid's net charge (0 → −1) never crosses zero → pI null
    assert summary.isoelectric_point is None
    assert any("isoelectric" in w for w in summary.warnings)


def test_amphoteric_pi_between_transitions() -> None:
    # +1 → 0 leg at pH 3, 0 → −1 leg at pH 9; pI between the two transitions.
    states = [
        Microstate(smiles="pos", charge=1, dg=dg_for_transition(3.0, 1)),
        Microstate(smiles="neut", charge=0, dg=0.0),
        Microstate(smiles="neg", charge=-1, dg=dg_for_transition(9.0, -1)),
    ]
    summary = compute_protonation_summary(
        states, working_ph=7.0, selected_forms=[("neut", 0.99), ("pos", 0.01)], ph_step=0.02
    )
    assert summary.isoelectric_point is not None
    assert 3.0 < summary.isoelectric_point < 9.0
    assert summary.pka_nearest_transition is not None
    assert min(abs(summary.pka_nearest_transition - t) for t in (3.0, 9.0)) < 0.05
    assert set(summary.charge_population) <= {-1, 0, 1}


def test_base_never_negative_isoelectric_null() -> None:
    states = [
        Microstate(smiles="B", charge=1, dg=0.0),
        Microstate(smiles="BH", charge=0, dg=5.0),
    ]
    summary = compute_protonation_summary(
        states, working_ph=7.0, selected_forms=[("B", 0.9), ("BH", 0.1)]
    )
    assert summary.isoelectric_point is None
    assert any("isoelectric" in w for w in summary.warnings)
    assert summary.charge_population[1] + summary.charge_population[0] == pytest.approx(1.0)


def test_empty_ensemble_returns_null_summary() -> None:
    summary = compute_protonation_summary([], working_ph=7.0, selected_forms=[])
    assert summary.top_two_occupancy_gap is None
    assert summary.occupancy_entropy is None
    assert summary.microstate_count == 0
    assert summary.warnings


def test_occupancies_reweight_uniform_states() -> None:
    states = [
        Microstate(smiles="X", charge=0, dg=0.0),
        Microstate(smiles="Y", charge=0, dg=0.0),
    ]
    assert compute_occupancies(states, ph=3.0) == pytest.approx([0.5, 0.5])


def test_charge_population_sums_to_one_for_multivalent() -> None:
    states = [
        Microstate(smiles="m1", charge=-2, dg=-1.0),
        Microstate(smiles="m2", charge=-1, dg=0.0),
        Microstate(smiles="m3", charge=-1, dg=1.0),
        Microstate(smiles="m4", charge=0, dg=2.0),
    ]
    occupancies = compute_occupancies(states, ph=7.0)
    assert sum(occupancies) == pytest.approx(1.0)
    summary = compute_protonation_summary(states, working_ph=7.0, selected_forms=[])
    assert sum(summary.charge_population.values()) == pytest.approx(1.0)
