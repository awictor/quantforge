"""Breeden-Litzenberger risk-neutral density implied by an SVI slice (svi_density)."""

import math

import pytest

from quantforge import svi_density
from quantforge.svi import SVIParams


S0, T, R = 100.0, 1.0, 0.02


def _flat_slice(sigma):
    # b=0 => w(k) = a = sigma^2 * t everywhere (flat smile).
    return SVIParams(a=sigma * sigma * T, b=0.0, rho=0.0, m=0.0, s=0.1)


def test_flat_slice_matches_lognormal():
    sigma = 0.2
    p = _flat_slice(sigma)
    F = S0 * math.exp(R * T)
    # Closed-form risk-neutral lognormal density of S_T at K.
    for K in (80.0, 100.0, 120.0):
        mu = math.log(F) - 0.5 * sigma * sigma * T
        sd = sigma * math.sqrt(T)
        pdf = math.exp(-((math.log(K) - mu) ** 2) / (2 * sd * sd)) / (
            K * sd * math.sqrt(2 * math.pi))
        g = svi_density(p, S0, T, R, K)
        assert g == pytest.approx(pdf, rel=2e-3)


def test_nonnegative_on_arbfree_slice():
    p = SVIParams(a=0.03, b=0.2, rho=-0.4, m=0.0, s=0.2)
    assert p.is_arbitrage_free_wings()
    for K in (60.0, 80.0, 100.0, 120.0, 150.0):
        assert svi_density(p, S0, T, R, K) >= 0.0


def test_integrates_to_one():
    p = SVIParams(a=0.03, b=0.15, rho=-0.3, m=0.0, s=0.25)
    # Trapezoid over a wide strike grid; RN density integrates to 1.
    lo, hi, n = 1.0, 500.0, 4000
    dK = (hi - lo) / n
    total = 0.0
    prev = svi_density(p, S0, T, R, lo)
    for i in range(1, n + 1):
        K = lo + i * dK
        cur = svi_density(p, S0, T, R, K)
        total += 0.5 * (prev + cur) * dK
        prev = cur
    assert total == pytest.approx(1.0, abs=5e-3)


def test_mean_equals_forward():
    # E[S_T] under the RN density = forward F.
    sigma = 0.25
    p = _flat_slice(sigma)
    F = S0 * math.exp(R * T)
    lo, hi, n = 1.0, 600.0, 6000
    dK = (hi - lo) / n
    mean = 0.0
    prev_K = lo
    prev = svi_density(p, S0, T, R, lo) * lo
    for i in range(1, n + 1):
        K = lo + i * dK
        cur = svi_density(p, S0, T, R, K) * K
        mean += 0.5 * (prev + cur) * dK
        prev = cur
    assert mean == pytest.approx(F, rel=3e-3)


def test_arbitrage_slice_can_go_negative():
    # Wing-bound-ok but butterfly-violating slice: g(k) < 0 near k ~ -0.6, so the
    # Breeden-Litzenberger density must be negative there.
    from quantforge.svi import svi_is_butterfly_free, svi_butterfly_arbitrage
    p = SVIParams(a=0.005, b=0.3, rho=-0.95, m=0.0, s=0.01)
    assert p.is_arbitrage_free_wings()
    assert not svi_is_butterfly_free(p)
    F = S0 * math.exp(R * T)
    bad_ks = svi_butterfly_arbitrage(p)
    vals = [svi_density(p, S0, T, R, F * math.exp(k)) for k in bad_ks]
    assert min(vals) < 0.0
