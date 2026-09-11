"""Cross-model parity: an SSVI surface fitted to SABR smiles reproduces the density.

Two independent risk-neutral density constructions -- the SABR
Breeden-Litzenberger density at an expiry and the SSVI-surface one calibrated to
SABR smiles at several expiries -- must agree where the fit is good.
"""

import math

import pytest

from quantforge.sabr import sabr_vol, sabr_density
from quantforge.ssvi import (
    calibrate_ssvi, ssvi_density, ssvi_is_arbitrage_free,
)


R = 0.02
ALPHA, BETA, RHO, NU = 0.2, 0.7, -0.3, 0.4
S0 = 100.0
EXPIRIES = [0.5, 1.0]
KS = [80.0, 90.0, 95.0, 100.0, 105.0, 110.0, 125.0]


def _fit_ssvi():
    market = []
    for t in EXPIRIES:
        F = S0 * math.exp(R * t)
        for K in KS:
            k = math.log(K / F)
            market.append((t, k, sabr_vol(F, K, t, ALPHA, BETA, RHO, NU)))
    return calibrate_ssvi(market)


def _sabr_dens(K, t):
    F = S0 * math.exp(R * t)
    return sabr_density(F, K, t, ALPHA, BETA, RHO, NU, r=R)


def test_ssvi_fits_sabr_surface_tightly():
    _p, rmse = _fit_ssvi()
    assert rmse < 5e-3


def test_surface_is_arbitrage_free():
    p, _ = _fit_ssvi()
    assert ssvi_is_arbitrage_free(p)


def test_densities_agree_near_the_money():
    p, _ = _fit_ssvi()
    t = 1.0
    for K in (95.0, 100.0, 105.0):
        assert ssvi_density(p, t, S0, R, K) == pytest.approx(
            _sabr_dens(K, t), rel=0.04)


def test_densities_agree_in_wings():
    p, _ = _fit_ssvi()
    t = 1.0
    for K in (90.0, 110.0):
        assert ssvi_density(p, t, S0, R, K) == pytest.approx(
            _sabr_dens(K, t), rel=0.06)


def test_both_densities_integrate_to_one():
    p, _ = _fit_ssvi()
    t = 1.0
    lo, hi, n = 1.0, 400.0, 4000
    dK = (hi - lo) / n
    tot_s = tot_v = 0.0
    ps = _sabr_dens(lo, t)
    pv = ssvi_density(p, t, S0, R, lo)
    for i in range(1, n + 1):
        K = lo + i * dK
        cs = _sabr_dens(K, t)
        cv = ssvi_density(p, t, S0, R, K)
        tot_s += 0.5 * (ps + cs) * dK
        tot_v += 0.5 * (pv + cv) * dK
        ps, pv = cs, cv
    assert tot_s == pytest.approx(1.0, abs=1e-2)
    assert tot_v == pytest.approx(1.0, abs=1e-2)
