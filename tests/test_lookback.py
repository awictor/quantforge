"""Tests for continuously-monitored lookback options.

Lookback closed forms are cross-checked against Monte Carlo out of band (slow),
so here we pin regression reference values and enforce the exact structural
identities the prices must satisfy:
  * expiry limits equal the realized intrinsic against the running extreme,
  * a lookback is worth at least the corresponding vanilla,
  * a floating-strike call/put is always in the money (payoff >= 0 by design).
"""

import math

import pytest

from quantforge import (
    floating_strike_lookback as fsl, fixed_strike_lookback as fxl,
    call_price, put_price, OptionType,
)


# --- Regression reference values (S=100, t=1, r=b=0.05, sigma=0.3) ---
def test_reference_values():
    assert fsl(100, 1.0, 0.05, 0.3, OptionType.CALL, b=0.05) == pytest.approx(23.788437, abs=1e-4)
    assert fsl(100, 1.0, 0.05, 0.3, OptionType.PUT, b=0.05) == pytest.approx(23.300731, abs=1e-4)
    assert fxl(100, 100, 1.0, 0.05, 0.3, OptionType.CALL, b=0.05) == pytest.approx(28.177788, abs=1e-4)
    assert fxl(100, 100, 1.0, 0.05, 0.3, OptionType.PUT, b=0.05) == pytest.approx(18.911379, abs=1e-4)


# --- Expiry limits ---
def test_floating_call_expiry_is_spot_minus_min():
    # At t->0 the floating call pays S - S_min.
    v = fsl(120, 1e-9, 0.05, 0.3, OptionType.CALL, s_extreme=100, b=0.05)
    assert v == pytest.approx(20.0, abs=1e-3)


def test_floating_put_expiry_is_max_minus_spot():
    v = fsl(100, 1e-9, 0.05, 0.3, OptionType.PUT, s_extreme=120, b=0.05)
    assert v == pytest.approx(20.0, abs=1e-3)


def test_fixed_call_expiry_is_max_minus_strike():
    v = fxl(130, 100, 1e-9, 0.05, 0.3, OptionType.CALL, s_extreme=130, b=0.05)
    assert v == pytest.approx(30.0, abs=1e-3)


def test_fixed_put_expiry_is_strike_minus_min():
    v = fxl(70, 100, 1e-9, 0.05, 0.3, OptionType.PUT, s_extreme=70, b=0.05)
    assert v == pytest.approx(30.0, abs=1e-3)


def test_fixed_call_otm_expiry_is_zero():
    v = fxl(100, 120, 1e-9, 0.05, 0.3, OptionType.CALL, s_extreme=100, b=0.05)
    assert v == pytest.approx(0.0, abs=1e-6)


# --- Lookback dominates vanilla ---
def test_lookback_call_at_least_vanilla():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.3
    van = call_price(S, K, t, r, sigma)
    assert fsl(S, t, r, sigma, OptionType.CALL, b=r) >= van - 1e-9
    assert fxl(S, K, t, r, sigma, OptionType.CALL, b=r) >= van - 1e-9


def test_lookback_put_at_least_vanilla():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.3
    van = put_price(S, K, t, r, sigma)
    assert fsl(S, t, r, sigma, OptionType.PUT, b=r) >= van - 1e-9
    assert fxl(S, K, t, r, sigma, OptionType.PUT, b=r) >= van - 1e-9


# --- Monotonicity ---
def test_fixed_call_decreasing_in_strike():
    kw = dict(S=100, t=1.0, r=0.05, sigma=0.3, option_type=OptionType.CALL, b=0.05)
    lo = fxl(K=90, **kw)
    hi = fxl(K=110, **kw)
    assert lo > hi


def test_prices_positive_across_vols():
    for sigma in (0.1, 0.25, 0.5, 0.8):
        assert fsl(100, 1.0, 0.05, sigma, OptionType.CALL, b=0.05) > 0
        assert fxl(100, 100, 1.0, 0.05, sigma, OptionType.PUT, b=0.05) > 0


# --- b == 0 (futures) is handled via the removable-singularity nudge ---
def test_zero_carry_is_finite():
    v = fsl(100, 1.0, 0.05, 0.3, OptionType.CALL, b=0.0)
    assert math.isfinite(v) and v > 0
