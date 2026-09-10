"""Tests for the Ho-Lee short-rate model."""

import math

import pytest

from quantforge import holee_zero_coupon_bond as zcb, holee_zero_coupon_yield as zcy


def test_bond_at_zero_is_one():
    assert zcb(0.03, 0.0, 0.002, 0.01) == 1.0


def test_short_yield_is_short_rate():
    assert zcy(0.03, 0.001, 0.002, 0.01) == pytest.approx(0.03, abs=1e-4)


def test_yield_matches_bond_price():
    r0, theta, sigma = 0.03, 0.002, 0.01
    for t in (1.0, 5.0, 10.0):
        assert zcy(r0, t, theta, sigma) == pytest.approx(
            -math.log(zcb(r0, t, theta, sigma)) / t, abs=1e-12)


def test_matches_monte_carlo():
    # Cross-checked against a Ho-Lee short-rate Monte Carlo (P(5) ~ 0.841).
    assert zcb(0.03, 5.0, 0.002, 0.01) == pytest.approx(0.841, abs=5e-3)


def test_convexity_lifts_bond():
    # Positive vol adds a +sigma^2 t^3/6 term, so the bond is worth more than
    # the zero-vol case.
    novol = zcb(0.03, 5.0, 0.002, 0.0)
    withvol = zcb(0.03, 5.0, 0.002, 0.05)
    assert withvol > novol


def test_higher_drift_lowers_bond():
    assert zcb(0.03, 5.0, 0.01, 0.01) < zcb(0.03, 5.0, 0.0, 0.01)


def test_rejects_negative_sigma():
    with pytest.raises(ValueError):
        zcb(0.03, 1.0, 0.002, -0.01)
