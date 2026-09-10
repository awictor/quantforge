"""Tests for the SVI butterfly (Gatheral-Jacquier g-function) arbitrage check."""

import math

import pytest

from quantforge import (
    SVIParams, svi_g, svi_butterfly_arbitrage, svi_is_butterfly_free,
    call_price, risk_neutral_density,
)


def test_reasonable_smile_is_butterfly_free():
    p = SVIParams(a=0.04, b=0.1, rho=-0.3, m=0.0, s=0.2)
    assert svi_is_butterfly_free(p)
    assert svi_butterfly_arbitrage(p) == []


def test_pathological_smile_flags_arbitrage():
    # Tiny level + steep, sharply-curved wings -> negative density region.
    p = SVIParams(a=0.001, b=0.9, rho=-0.9, m=0.0, s=0.02)
    assert not svi_is_butterfly_free(p)
    assert len(svi_butterfly_arbitrage(p)) > 0


def test_g_matches_breeden_litzenberger_density_sign():
    # Where g >= 0 the extracted risk-neutral density must be non-negative.
    p = SVIParams(a=0.04, b=0.4, rho=-0.6, m=0.0, s=0.1)
    t, F = 1.0, 100.0
    strikes = [k * 0.5 for k in range(2, 600)]

    def vol(K):
        w = p.total_variance(math.log(K / F))
        return math.sqrt(max(w, 1e-9) / t)

    calls = [call_price(F, K, t, 0.0, vol(K)) for K in strikes]
    _, dens = risk_neutral_density(strikes, calls, t, 0.0)
    assert all(d >= -1e-9 for d in dens)
    assert all(svi_g(p, math.log(K / F)) >= -1e-9 for K in strikes)


def test_g_positive_at_the_money_for_typical_smile():
    p = SVIParams(a=0.04, b=0.2, rho=-0.4, m=0.0, s=0.3)
    assert svi_g(p, 0.0) > 0


def test_flat_smile_is_butterfly_free():
    # b = 0: flat total variance -> lognormal -> always arbitrage free.
    p = SVIParams(a=0.04, b=0.0, rho=0.0, m=0.0, s=0.1)
    assert svi_is_butterfly_free(p)
