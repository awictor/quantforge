"""Cross-model parity: an SVI slice fitted to a SABR smile reproduces its density.

Two independent risk-neutral density constructions -- the SABR
Breeden-Litzenberger density and the SVI one on a slice calibrated to the same
smile -- must agree where the fit is good.
"""

import math

import pytest

from quantforge.sabr import sabr_vol, sabr_density
from quantforge.svi import calibrate_svi, svi_density


F, T, R = 100.0, 1.0, 0.02
ALPHA, BETA, RHO, NU = 0.2, 0.7, -0.3, 0.4
S0 = F * math.exp(-R * T)
KS = [70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 135.0]


def _fit_svi():
    ks = [math.log(K / F) for K in KS]
    tv = [sabr_vol(F, K, T, ALPHA, BETA, RHO, NU) ** 2 * T for K in KS]
    return calibrate_svi(ks, tv)


def test_svi_fits_sabr_smile_tightly():
    _p, rmse = _fit_svi()
    assert rmse < 1e-3


def test_densities_agree_near_the_money():
    p, _ = _fit_svi()
    for K in (95.0, 100.0, 105.0):
        ds = sabr_density(F, K, T, ALPHA, BETA, RHO, NU, r=R)
        dv = svi_density(p, S0, T, R, K)
        assert dv == pytest.approx(ds, rel=0.02)


def test_densities_agree_in_wings():
    p, _ = _fit_svi()
    for K in (85.0, 115.0):
        ds = sabr_density(F, K, T, ALPHA, BETA, RHO, NU, r=R)
        dv = svi_density(p, S0, T, R, K)
        assert dv == pytest.approx(ds, rel=0.05)


def test_both_densities_nonnegative():
    p, _ = _fit_svi()
    for K in KS:
        assert sabr_density(F, K, T, ALPHA, BETA, RHO, NU, r=R) >= 0.0
        assert svi_density(p, S0, T, R, K) >= 0.0


def test_both_integrate_to_one():
    p, _ = _fit_svi()
    lo, hi, n = 1.0, 400.0, 4000
    dK = (hi - lo) / n
    tot_s = tot_v = 0.0
    ps = sabr_density(F, lo, T, ALPHA, BETA, RHO, NU, r=R)
    pv = svi_density(p, S0, T, R, lo)
    for i in range(1, n + 1):
        K = lo + i * dK
        cs = sabr_density(F, K, T, ALPHA, BETA, RHO, NU, r=R)
        cv = svi_density(p, S0, T, R, K)
        tot_s += 0.5 * (ps + cs) * dK
        tot_v += 0.5 * (pv + cv) * dK
        ps, pv = cs, cv
    assert tot_s == pytest.approx(1.0, abs=1e-2)
    assert tot_v == pytest.approx(1.0, abs=1e-2)
