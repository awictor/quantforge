"""Coupon-bond analytics: price, YTM, duration, convexity, DV01."""

import pytest

from quantforge import (
    bond_cashflows, bond_price_from_yield, macaulay_duration,
    modified_duration, convexity, bond_dv01, yield_to_maturity,
)


CF = bond_cashflows(100.0, 0.05, 5.0, freq=2)
Y = 0.04


def test_cashflow_schedule():
    assert len(CF) == 10
    assert CF[-1] == (5.0, 102.5)   # last coupon + face
    assert CF[0] == (0.5, 2.5)


def test_ytm_round_trip():
    p = bond_price_from_yield(CF, Y)
    assert yield_to_maturity(CF, p) == pytest.approx(Y, abs=1e-8)


def test_duration_matches_finite_difference():
    p = bond_price_from_yield(CF, Y)
    h = 1e-6
    fd = -(bond_price_from_yield(CF, Y + h) - bond_price_from_yield(CF, Y - h)) / (2 * h) / p
    assert modified_duration(CF, Y) == pytest.approx(fd, abs=1e-4)


def test_convexity_matches_finite_difference():
    p = bond_price_from_yield(CF, Y)
    h = 1e-6
    fd = (bond_price_from_yield(CF, Y + h) - 2 * p + bond_price_from_yield(CF, Y - h)) / (h * h) / p
    assert convexity(CF, Y) == pytest.approx(fd, rel=1e-3)


def test_duration_below_maturity_and_convexity_positive():
    assert 0 < macaulay_duration(CF, Y) < 5.0
    assert convexity(CF, Y) > 0.0


def test_continuous_modified_equals_macaulay():
    assert modified_duration(CF, Y) == macaulay_duration(CF, Y)


def test_dv01_sign_and_magnitude():
    p = bond_price_from_yield(CF, Y)
    dv = bond_dv01(CF, Y)
    assert dv < 0.0
    assert dv == pytest.approx(-modified_duration(CF, Y) * p * 1e-4, abs=1e-12)


def test_zero_coupon_duration_is_maturity():
    zc = [(5.0, 100.0)]
    assert macaulay_duration(zc, Y) == pytest.approx(5.0, abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        yield_to_maturity(CF, -1.0)
    with pytest.raises(ValueError):
        bond_cashflows(100.0, 0.05, -1.0)
    with pytest.raises(ValueError):
        yield_to_maturity(CF, 1e9)  # price above sum of cashflows
