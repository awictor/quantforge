"""Surface-consistent smile consumers on an SSVI slice.

Variance-swap strike, VIX/SVIX, Breeden-Litzenberger density, and BKM moments
built from the surface's smile at a fitted expiry.
"""

import math

import pytest

from quantforge import (
    ssvi_variance_swap_strike, ssvi_vix, ssvi_svix, ssvi_density,
    ssvi_bkm_moments,
)
from quantforge.ssvi import SSVIParams, ssvi_is_arbitrage_free


S0, R, SIGMA = 100.0, 0.02, 0.2
TS = [0.25, 0.5, 1.0]


def _surface(rho=-0.3, eta=0.5, gamma=0.4):
    # theta_t = sigma^2 * t keeps the ATM vol flat across expiries.
    return SSVIParams(rho=rho, eta=eta, gamma=gamma,
                      thetas={t: SIGMA * SIGMA * t for t in TS})


def test_surface_is_arbitrage_free():
    assert ssvi_is_arbitrage_free(_surface())


def test_varswap_at_least_atm_variance():
    # With a skew (eta>0) the fair variance exceeds the flat ATM variance.
    vs = ssvi_variance_swap_strike(_surface(), 1.0, S0, R)
    assert vs > SIGMA * SIGMA - 1e-4
    assert math.sqrt(vs) == pytest.approx(SIGMA, abs=0.02)


def test_vix_near_atm_vol():
    vix = ssvi_vix(_surface(), 1.0, S0, R)
    assert vix == pytest.approx(100 * SIGMA, abs=1.5)


def test_svix_positive_and_near_vol():
    svix = ssvi_svix(_surface(), 1.0, S0, R)
    assert svix > 0.0
    assert svix == pytest.approx(100 * SIGMA, abs=2.0)


def test_density_integrates_to_one_and_mean_is_forward():
    p = _surface()
    t = 1.0
    F = S0 * math.exp(R * t)
    lo, hi, n = 1.0, 500.0, 4000
    dK = (hi - lo) / n
    total = mean = 0.0
    prev = ssvi_density(p, t, S0, R, lo)
    prev_m = prev * lo
    for i in range(1, n + 1):
        K = lo + i * dK
        cur = ssvi_density(p, t, S0, R, K)
        total += 0.5 * (prev + cur) * dK
        mean += 0.5 * (prev_m + cur * K) * dK
        prev, prev_m = cur, cur * K
    assert total == pytest.approx(1.0, abs=5e-3)
    assert mean == pytest.approx(F, rel=3e-3)


def test_density_nonnegative_on_arbfree_slice():
    p = _surface()
    F = S0 * math.exp(R)
    for k in (-0.5, -0.2, 0.0, 0.2, 0.5):
        assert ssvi_density(p, 1.0, S0, R, F * math.exp(k)) >= 0.0


def test_bkm_skew_sign_follows_rho():
    _v, skew_neg, kurt = ssvi_bkm_moments(_surface(rho=-0.4), 1.0, S0, R)
    assert skew_neg < 0.0
    assert kurt > 0.0
    _v2, skew_pos, _k = ssvi_bkm_moments(_surface(rho=0.4), 1.0, S0, R)
    assert skew_pos > 0.0


def test_unfitted_expiry_raises():
    with pytest.raises(ValueError):
        ssvi_density(_surface(), 0.75, S0, R, 100.0)
