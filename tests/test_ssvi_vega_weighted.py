"""Vega-weighted SSVI calibration favours the near-the-money fit."""

import math
import random

import pytest

from quantforge import SSVIParams, calibrate_ssvi, ssvi_is_arbitrage_free


TRUE = SSVIParams(rho=-0.4, eta=1.0, gamma=0.5,
                  thetas={0.25: 0.01, 1.0: 0.045, 2.0: 0.10})
KS = [-0.4, -0.2, 0.0, 0.2, 0.4]


def _clean_market():
    return [(t, k, math.sqrt(TRUE.total_variance(k, t) / t))
            for t in TRUE.thetas for k in KS]


def _wing_noisy_market(seed=1):
    rng = random.Random(seed)
    m = []
    for t in TRUE.thetas:
        for k in KS:
            iv = math.sqrt(TRUE.total_variance(k, t) / t)
            noise = 0.01 * abs(k) / 0.4    # larger noise in the wings
            m.append((t, k, iv + rng.uniform(-noise, noise)))
    return m


def _atm_error(params):
    e = 0.0
    for t in TRUE.thetas:
        got = params.implied_vol(0.0, t)
        ref = math.sqrt(TRUE.total_variance(0.0, t) / t)
        e += abs(got - ref)
    return e


def test_clean_surface_recovered_with_vega_weighting():
    params, rmse = calibrate_ssvi(_clean_market(), vega_weighted=True)
    assert rmse < 1e-4
    assert params.rho == pytest.approx(TRUE.rho, abs=1e-2)


def test_vega_weighting_improves_atm_fit_under_wing_noise():
    market = _wing_noisy_market()
    p_w, _ = calibrate_ssvi(market, vega_weighted=True)
    p_u, _ = calibrate_ssvi(market, vega_weighted=False)
    assert _atm_error(p_w) < _atm_error(p_u)


def test_vega_weighted_stays_valid():
    params, _ = calibrate_ssvi(_clean_market(), vega_weighted=True)
    assert params.eta > 0 and 0 < params.gamma < 1 and -1 < params.rho < 1
    assert ssvi_is_arbitrage_free(params)


def test_default_is_unweighted():
    # Without the flag the fit is the plain (unweighted) least squares.
    market = _clean_market()
    p_default, r_default = calibrate_ssvi(market)
    p_unw, r_unw = calibrate_ssvi(market, vega_weighted=False)
    assert r_default == pytest.approx(r_unw, abs=1e-9)
