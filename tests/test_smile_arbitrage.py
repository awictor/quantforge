"""Model-free butterfly-arbitrage check on a smile (smile_arbitrage_violations)."""

import math

import pytest

from quantforge import smile_arbitrage_violations, smile_is_arbitrage_free
from quantforge.svi import SVIParams, svi_is_butterfly_free


S0, R, T = 100.0, 0.02, 1.0
F = S0 * math.exp(R * T)


def _svi_vol(p):
    return lambda K: p.implied_vol(math.log(K / F), T)


def test_flat_smile_is_arbitrage_free():
    assert smile_is_arbitrage_free(S0, T, R, lambda K: 0.2)
    assert smile_arbitrage_violations(S0, T, R, lambda K: 0.2) == []


def test_convex_svi_slice_is_arbitrage_free():
    p = SVIParams(a=0.03, b=0.2, rho=-0.3, m=0.0, s=0.25)
    assert svi_is_butterfly_free(p)                       # SVI's own g-function
    assert smile_is_arbitrage_free(S0, T, R, _svi_vol(p))  # density-based check


def test_butterfly_violating_slice_is_flagged():
    # Steep, high-curvature slice with g(k) < 0 near k ~ -0.6 (see loop 305).
    p = SVIParams(a=0.005, b=0.3, rho=-0.95, m=0.0, s=0.01)
    assert not svi_is_butterfly_free(p)
    assert not smile_is_arbitrage_free(S0, T, R, _svi_vol(p))
    assert len(smile_arbitrage_violations(S0, T, R, _svi_vol(p))) > 0


def test_agrees_with_svi_g_function_across_slices():
    # The density-sign check and SVI's analytic g-function should agree on
    # whether each slice is butterfly-free.
    slices = [
        SVIParams(a=0.03, b=0.2, rho=-0.3, m=0.0, s=0.25),   # free
        SVIParams(a=0.04, b=0.15, rho=-0.2, m=0.0, s=0.2),   # free
        SVIParams(a=0.005, b=0.3, rho=-0.95, m=0.0, s=0.01),  # violating
    ]
    for p in slices:
        assert smile_is_arbitrage_free(S0, T, R, _svi_vol(p)) == \
            svi_is_butterfly_free(p)


def test_flat_smile_with_dividends():
    # A flat smile stays arbitrage-free with a nonzero dividend yield.
    assert smile_is_arbitrage_free(S0, T, R, lambda K: 0.25, q=0.03)
