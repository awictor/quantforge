"""Tests for the Cox-Ingersoll-Ross short-rate model."""

import math

import pytest

from quantforge import cir_zero_coupon_bond as zcb, cir_zero_coupon_yield as zcy


KTS = dict(kappa=0.5, theta=0.05, sigma=0.08)


def test_bond_at_zero_is_one():
    assert zcb(0.03, 0.0, **KTS) == 1.0


def test_short_yield_is_short_rate():
    assert zcy(0.03, 0.01, **KTS) == pytest.approx(0.03, abs=1e-3)


def test_long_yield_near_cir_limit():
    g = math.sqrt(0.5 ** 2 + 2 * 0.08 ** 2)
    limit = 2 * 0.5 * 0.05 / (g + 0.5)
    assert zcy(0.03, 50.0, **KTS) == pytest.approx(limit, abs=2e-3)


def test_bond_decreasing_in_maturity():
    assert zcb(0.03, 5.0, **KTS) < zcb(0.03, 1.0, **KTS) < 1.0


def test_matches_monte_carlo():
    # Cross-checked against a floored-Euler CIR Monte Carlo (P(5) ~ 0.809).
    assert zcb(0.03, 5.0, **KTS) == pytest.approx(0.809, abs=5e-3)


def test_higher_rate_lowers_bond():
    assert zcb(0.06, 5.0, **KTS) < zcb(0.02, 5.0, **KTS)


def test_rejects_negative_rate():
    with pytest.raises(ValueError):
        zcb(-0.01, 1.0, **KTS)


def test_rejects_bad_params():
    with pytest.raises(ValueError):
        zcb(0.03, 1.0, kappa=0.0, theta=0.05, sigma=0.08)
