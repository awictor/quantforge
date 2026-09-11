"""Variance-swap strike replicated from an SVI slice (svi_variance_swap_strike)."""

import pytest

from quantforge import svi_variance_swap_strike, SVIParams


S0, T, R = 100.0, 1.0, 0.03


def test_flat_slice_returns_flat_variance():
    # b = 0 -> flat total variance a = sigma^2 t -> flat vol sigma; the swap
    # strike is that variance (up to strip truncation).
    sig = 0.2
    flat = SVIParams(a=sig * sig * T, b=0.0, rho=0.0, m=0.0, s=0.1)
    strike = svi_variance_swap_strike(flat, S0, T, R, n_strikes=801, width=10.0)
    assert strike == pytest.approx(sig * sig, abs=2e-3)


def test_skew_lifts_strike_above_atm_variance():
    # A skewed/convex slice makes the fair variance exceed the ATM variance.
    p = SVIParams(a=0.03, b=0.1, rho=-0.5, m=0.0, s=0.2)
    atm_var = p.implied_vol(0.0, T) ** 2
    strike = svi_variance_swap_strike(p, S0, T, R)
    assert strike > atm_var


def test_positive():
    p = SVIParams(a=0.03, b=0.1, rho=-0.5, m=0.0, s=0.2)
    assert svi_variance_swap_strike(p, S0, T, R) > 0.0


def test_dividend_carry_runs():
    p = SVIParams(a=0.03, b=0.08, rho=-0.3, m=0.0, s=0.25)
    v = svi_variance_swap_strike(p, S0, T, R, q=0.02)
    assert v > 0.0


def test_higher_convexity_raises_strike():
    lo = SVIParams(a=0.04, b=0.05, rho=0.0, m=0.0, s=0.3)
    hi = SVIParams(a=0.04, b=0.15, rho=0.0, m=0.0, s=0.3)
    assert svi_variance_swap_strike(hi, S0, T, R) > svi_variance_swap_strike(lo, S0, T, R)
