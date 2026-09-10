"""Tests for composite (compo) FX options."""

import math

import pytest

from quantforge import compo_option, call_price, OptionType


BASE = dict(S=100, K=100, t=1.0, r_domestic=0.05, r_foreign=0.05,
            sigma_asset=0.2, sigma_fx=0.1)


def test_zero_fx_vol_is_black_scholes():
    # No FX vol -> the combined vol is just the asset vol -> plain BSM.
    v = compo_option(**{**BASE, "sigma_fx": 0.0}, rho=0.0, option_type=OptionType.CALL)
    assert v == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2), abs=1e-9)


def test_higher_correlation_raises_price():
    lo = compo_option(**BASE, rho=-0.5, option_type=OptionType.CALL)
    mid = compo_option(**BASE, rho=0.0, option_type=OptionType.CALL)
    hi = compo_option(**BASE, rho=0.8, option_type=OptionType.CALL)
    assert lo < mid < hi


def test_matches_monte_carlo():
    # Cross-checked against a combined-lognormal Monte Carlo (~12.89 at rho=0.5).
    v = compo_option(**BASE, rho=0.5, option_type=OptionType.CALL)
    assert v == pytest.approx(12.89, abs=0.1)


def test_combined_vol_formula():
    # The effective vol is sqrt(sa^2 + sfx^2 + 2 rho sa sfx); at rho=1 it is the
    # sum sa + sfx, so the compo call equals a BSM call at that vol.
    v = compo_option(**BASE, rho=1.0, option_type=OptionType.CALL)
    assert v == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2 + 0.1), abs=1e-9)


def test_put_call_parity():
    c = compo_option(**BASE, rho=0.3, option_type=OptionType.CALL)
    p = compo_option(**BASE, rho=0.3, option_type=OptionType.PUT)
    assert c - p == pytest.approx(100 - 100 * math.exp(-0.05), abs=1e-9)


def test_rejects_bad_correlation():
    with pytest.raises(ValueError):
        compo_option(**BASE, rho=1.5)
