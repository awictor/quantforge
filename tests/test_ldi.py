"""Liability-driven investing: funding, duration hedging, surplus risk."""

import math

import pytest

from quantforge import (
    liability_pv, funding_ratio, surplus, liability_duration, hedge_ratio,
    required_hedge_duration, surplus_at_risk,
    liability_convexity, surplus_change_under_shock, funded_ratio_return,
)
from quantforge.mathfns import norm_ppf


CF = [(5, 100), (10, 150), (15, 200), (20, 250), (30, 300)]
R = 0.03


def test_funding_ratio_and_surplus():
    L = liability_pv(CF, R)
    assert funding_ratio(L * 1.1, L) == pytest.approx(1.1)
    assert surplus(L * 1.1, L) == pytest.approx(0.1 * L, abs=1e-6)


def test_liability_duration_long_dated():
    assert liability_duration(CF, R) > 10


def test_required_hedge_duration_gives_unit_ratio():
    L = liability_pv(CF, R)
    D = liability_duration(CF, R)
    A = L * 1.1
    req = required_hedge_duration(A, D, L)
    assert hedge_ratio(req, A, D, L) == pytest.approx(1.0)


def test_hedge_ratio_scales():
    L = liability_pv(CF, R)
    D = liability_duration(CF, R)
    A = L * 1.1
    req = required_hedge_duration(A, D, L)
    assert hedge_ratio(2 * req, A, D, L) == pytest.approx(2.0)


def test_surplus_at_risk_positive_and_monotone():
    L = liability_pv(CF, R)
    s95 = surplus_at_risk(L * 1.1, L, 0.08, 0.95)
    s99 = surplus_at_risk(L * 1.1, L, 0.08, 0.99)
    assert s95 > 0
    assert s99 > s95
    assert s95 == pytest.approx(norm_ppf(0.95) * 0.08 * L, abs=1e-6)


def test_convexity_positive():
    assert liability_convexity(CF, R) > 0


def test_immunized_surplus_stable():
    L = liability_pv(CF, R)
    D = liability_duration(CF, R)
    C = liability_convexity(CF, R)
    # Matched duration and convexity: surplus unchanged under a shock.
    assert surplus_change_under_shock(L, D, C, L, D, C, 0.01) == pytest.approx(0.0, abs=1e-9)


def test_duration_only_moves_at_second_order():
    L = liability_pv(CF, R)
    D = liability_duration(CF, R)
    C = liability_convexity(CF, R)
    matched = abs(surplus_change_under_shock(L, D, C, L, D, C, 0.01))
    dur_only = abs(surplus_change_under_shock(L, D, 0.5 * C, L, D, C, 0.01))
    assert dur_only > matched


def test_funded_ratio_return():
    assert funded_ratio_return(0.08, 0.05, 1.1) == pytest.approx(1.08 / 1.05 - 1)
    # Independent of the funding-ratio level.
    assert funded_ratio_return(0.08, 0.05, 1.1) == pytest.approx(
        funded_ratio_return(0.08, 0.05, 0.9))


def test_convexity_validation():
    with pytest.raises(ValueError):
        funded_ratio_return(0.08, -1.5, 1.1)


def test_validation():
    L = liability_pv(CF, R)
    with pytest.raises(ValueError):
        funding_ratio(L, 0)
    with pytest.raises(ValueError):
        surplus_at_risk(L, L, 0.08, 0.3)
