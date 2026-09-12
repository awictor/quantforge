"""Equity valuation: cost of capital, dividend discount, DCF."""

import pytest

from quantforge import (
    capm_cost_of_equity, wacc, gordon_growth_value, terminal_value, two_stage_dcf,
)


def test_capm():
    assert capm_cost_of_equity(0.03, 1.2, 0.05) == pytest.approx(0.03 + 1.2 * 0.05)


def test_wacc_between_costs_and_all_equity():
    ke, kd, tax = 0.10, 0.05, 0.21
    w = wacc(600, 400, ke, kd, tax)
    assert kd * (1 - tax) < w < ke
    assert wacc(1000, 0, ke, kd, tax) == pytest.approx(ke)


def test_gordon_growth():
    assert gordon_growth_value(2, 0.08, 0.03) == pytest.approx(40)


def test_terminal_value():
    assert terminal_value(100, 0.08, 0.03) == pytest.approx(100 * 1.03 / 0.05)


def test_two_stage_dcf_monotonicity():
    v = two_stage_dcf([100] * 5, 0.08, 0.03)
    assert v > 0
    assert two_stage_dcf([100] * 5, 0.08, 0.04) > v   # higher growth
    assert two_stage_dcf([100] * 5, 0.10, 0.03) < v   # higher discount


def test_validation():
    with pytest.raises(ValueError):
        gordon_growth_value(2, 0.03, 0.05)   # r <= g
    with pytest.raises(ValueError):
        two_stage_dcf([100], 0.03, 0.05)     # r <= g
