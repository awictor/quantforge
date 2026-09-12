"""Brinson performance attribution."""

import pytest

from quantforge import (
    allocation_effect, selection_effect, interaction_effect, brinson_attribution,
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


def test_validation():
    with pytest.raises(ValueError):
        allocation_effect([0.5], [0.4, 0.4], [0.08, 0.06])
    with pytest.raises(ValueError):
        selection_effect([], [], [])
