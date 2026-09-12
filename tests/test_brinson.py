"""Brinson performance attribution."""

import pytest

from quantforge import (
    allocation_effect, selection_effect, interaction_effect, brinson_attribution,
    carino_factor, linked_active_return, carino_linked_effects,
)


PW = [0.5, 0.3, 0.2]
BW = [0.4, 0.4, 0.2]
PR = [0.10, 0.05, 0.08]
BR = [0.08, 0.06, 0.07]


def test_effects_sum_to_active_return():
    res = brinson_attribution(PW, BW, PR, BR)
    assert res["total_effect"] == pytest.approx(res["active_return"], abs=1e-12)


def test_effect_formulas():
    res = brinson_attribution(PW, BW, PR, BR)
    assert res["allocation"][0] == pytest.approx((0.5 - 0.4) * 0.08)
    assert res["selection"][0] == pytest.approx(0.4 * (0.10 - 0.08))
    assert res["interaction"][0] == pytest.approx((0.5 - 0.4) * (0.10 - 0.08))


def test_identical_portfolio_zero():
    res = brinson_attribution(BW, BW, BR, BR)
    assert res["active_return"] == pytest.approx(0.0)
    assert res["total_effect"] == pytest.approx(0.0)


def test_pure_allocation_isolates_allocation():
    res = brinson_attribution(PW, BW, BR, BR)  # same returns, different weights
    assert sum(res["selection"]) == pytest.approx(0.0)
    assert sum(res["interaction"]) == pytest.approx(0.0)
    assert abs(sum(res["allocation"])) > 0


RP = [0.10, 0.05, -0.02]
RB = [0.08, 0.06, 0.00]


def test_carino_factor_equal_returns():
    assert carino_factor(0.05, 0.05) == pytest.approx(1 / 1.05)


def test_linked_active_return_geometric():
    pp = 1.10 * 1.05 * 0.98 - 1
    pb = 1.08 * 1.06 * 1.00 - 1
    assert linked_active_return(RP, RB) == pytest.approx(pp - pb)


def test_carino_linked_effects_sum_to_linked_active():
    per_active = [RP[i] - RB[i] for i in range(3)]
    smoothed = carino_linked_effects(per_active, RP, RB)
    assert sum(smoothed) == pytest.approx(linked_active_return(RP, RB), abs=1e-10)


def test_arithmetic_sum_leaves_residual():
    per_active = [RP[i] - RB[i] for i in range(3)]
    assert abs(sum(per_active) - linked_active_return(RP, RB)) > 1e-6


def test_single_period_linked_is_arithmetic():
    assert linked_active_return([0.05], [0.03]) == pytest.approx(0.05 - 0.03)


def test_carino_validation():
    with pytest.raises(ValueError):
        carino_factor(-1.5, 0.03)


def test_validation():
    with pytest.raises(ValueError):
        allocation_effect([0.5], [0.4, 0.4], [0.08, 0.06])
    with pytest.raises(ValueError):
        selection_effect([], [], [])
