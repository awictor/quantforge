"""Tests for the Vasicek short-rate model."""

import math

import pytest

from quantforge import zero_coupon_bond, zero_coupon_yield, bond_option, OptionType


KTS = dict(kappa=0.5, theta=0.05, sigma=0.01)


def test_bond_price_at_zero_is_one():
    assert zero_coupon_bond(0.03, 0.0, **KTS) == 1.0


def test_short_yield_is_short_rate():
    y = zero_coupon_yield(0.03, 0.01, **KTS)
    assert y == pytest.approx(0.03, abs=1e-3)


def test_long_yield_near_long_run_minus_convexity():
    y = zero_coupon_yield(0.03, 30.0, **KTS)
    limit = 0.05 - 0.01 ** 2 / (2 * 0.5 ** 2)
    assert y == pytest.approx(limit, abs=2e-3)


def test_bond_price_decreasing_in_maturity():
    p1 = zero_coupon_bond(0.03, 1.0, **KTS)
    p5 = zero_coupon_bond(0.03, 5.0, **KTS)
    assert p5 < p1 < 1.0


def test_bond_option_put_call_parity():
    r0, K = 0.03, 0.85
    c = bond_option(r0, 1.0, 5.0, K, **KTS, option_type=OptionType.CALL)
    p = bond_option(r0, 1.0, 5.0, K, **KTS, option_type=OptionType.PUT)
    Pb = zero_coupon_bond(r0, 5.0, **KTS)
    Po = zero_coupon_bond(r0, 1.0, **KTS)
    assert c - p == pytest.approx(Pb - K * Po, abs=1e-9)


def test_bond_option_matches_monte_carlo():
    # Cross-checked against a short-rate Monte Carlo (~0.0007).
    c = bond_option(0.03, 1.0, 5.0, 0.85, **KTS, option_type=OptionType.CALL)
    assert c == pytest.approx(0.0007, abs=5e-4)


def test_higher_vol_raises_bond_option():
    lo = bond_option(0.03, 1.0, 5.0, 0.88, kappa=0.5, theta=0.05, sigma=0.005)
    hi = bond_option(0.03, 1.0, 5.0, 0.88, kappa=0.5, theta=0.05, sigma=0.02)
    assert hi > lo


def test_rejects_bad_option_maturities():
    with pytest.raises(ValueError):
        bond_option(0.03, 5.0, 1.0, 0.85, **KTS)   # option after bond
