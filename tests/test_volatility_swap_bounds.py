"""Volatility-swap strike bounds from a smile (volatility_swap_bounds_from_smile)."""

import math

import pytest

from quantforge import (
    volatility_swap_bounds_from_smile, variance_swap_from_smile,
)
from quantforge.svi import SVIParams


S0, R, T = 100.0, 0.02, 1.0


def _svi_vol(p):
    F = S0 * math.exp(R * T)
    return lambda K: p.implied_vol(math.log(K / F), T)


def test_flat_smile_collapses_bounds():
    flat = lambda K: 0.2
    lo, hi = volatility_swap_bounds_from_smile(S0, T, R, flat)
    # Finite-strip replication of a flat smile slightly undershoots the exact
    # variance, so both bounds sit a touch below 0.20; the bracket still
    # collapses (no convexity premium for a flat smile).
    assert lo == pytest.approx(0.2, abs=5e-3)
    assert hi == pytest.approx(0.2, abs=5e-3)
    assert hi - lo == pytest.approx(0.0, abs=1e-3)


def test_convex_smile_orders_lower_below_upper():
    p = SVIParams(a=0.03, b=0.2, rho=-0.3, m=0.0, s=0.25)
    lo, hi = volatility_swap_bounds_from_smile(S0, T, R, _svi_vol(p))
    assert lo < hi                       # a convexity premium exists
    assert hi > 0.0 and lo > 0.0


def test_upper_is_sqrt_variance_strike():
    p = SVIParams(a=0.04, b=0.15, rho=-0.2, m=0.0, s=0.2)
    vf = _svi_vol(p)
    _lo, hi = volatility_swap_bounds_from_smile(S0, T, R, vf)
    var = variance_swap_from_smile(S0, T, R, vf)
    assert hi == pytest.approx(math.sqrt(var), rel=1e-9)


def test_lower_is_atmf_vol():
    p = SVIParams(a=0.03, b=0.2, rho=-0.4, m=0.0, s=0.2)
    vf = _svi_vol(p)
    F = S0 * math.exp(R * T)
    lo, _hi = volatility_swap_bounds_from_smile(S0, T, R, vf)
    assert lo == pytest.approx(vf(F), rel=1e-9)


def test_convexity_premium_grows_with_vol_of_vol():
    # A wider SVI curvature (larger b) -> more convex smile -> wider bracket.
    narrow = SVIParams(a=0.03, b=0.10, rho=-0.3, m=0.0, s=0.25)
    wide = SVIParams(a=0.03, b=0.30, rho=-0.3, m=0.0, s=0.25)
    lo_n, hi_n = volatility_swap_bounds_from_smile(S0, T, R, _svi_vol(narrow))
    lo_w, hi_w = volatility_swap_bounds_from_smile(S0, T, R, _svi_vol(wide))
    assert (hi_w - lo_w) > (hi_n - lo_n)
