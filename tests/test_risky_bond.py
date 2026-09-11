"""Defaultable (risky) coupon-bond pricing (credit.risky_bond_price)."""

import math

import pytest

from quantforge import (
    SurvivalCurve, risky_bond_price, risky_bond_yield_spread,
)
from quantforge.bondmath import bond_cashflows, bond_price_from_yield


CF = bond_cashflows(100.0, 0.05, 5.0, freq=2)
R = 0.03


def _curve(h):
    return SurvivalCurve([1, 3, 5, 10], [h] * 4)


def test_zero_hazard_equals_risk_free():
    rf = bond_price_from_yield(CF, R)
    assert risky_bond_price(_curve(0.0), CF, R, 0.4) == pytest.approx(rf, abs=1e-6)


def test_risky_below_risk_free():
    rf = bond_price_from_yield(CF, R)
    assert risky_bond_price(_curve(0.03), CF, R, 0.4) < rf


def test_higher_hazard_cheaper():
    lo = risky_bond_price(_curve(0.03), CF, R, 0.4)
    hi = risky_bond_price(_curve(0.06), CF, R, 0.4)
    assert hi < lo


def test_higher_recovery_richer():
    low_rec = risky_bond_price(_curve(0.04), CF, R, 0.2)
    high_rec = risky_bond_price(_curve(0.04), CF, R, 0.6)
    assert high_rec > low_rec


def test_yield_spread_positive_and_reprices():
    c = _curve(0.03)
    s = risky_bond_yield_spread(c, CF, R, 0.4)
    assert s > 0.0
    target = risky_bond_price(c, CF, R, 0.4)
    reprice = sum(cf * math.exp(-(R + s) * t) for t, cf in CF)
    assert reprice == pytest.approx(target, abs=1e-6)


def test_zero_hazard_zero_spread():
    s = risky_bond_yield_spread(_curve(0.0), CF, R, 0.4)
    assert s == pytest.approx(0.0, abs=1e-6)
