"""Arbitrage-free SSVI calibration via a ramped no-arb penalty."""

import math

import pytest

from quantforge import (
    SSVIParams,
    calibrate_ssvi,
    calibrate_ssvi_arbitrage_free,
    ssvi_is_arbitrage_free,
)


def _market_from(params):
    m = []
    for t, th in params.thetas.items():
        for k in (-0.3, -0.15, 0.0, 0.15, 0.3):
            w = params.total_variance(k, t)
            m.append((t, k, math.sqrt(w / t)))
    return m


ARB = SSVIParams(rho=-0.5, eta=8.0, gamma=0.3,
                 thetas={0.1: 0.02, 0.5: 0.03, 1.0: 0.05})
GOOD = SSVIParams(rho=-0.4, eta=1.0, gamma=0.5,
                  thetas={0.25: 0.01, 1.0: 0.045, 2.0: 0.10})


def test_constructed_arb_surface_is_flagged():
    assert not ssvi_is_arbitrage_free(ARB)


def test_arb_free_fit_produces_arbitrage_free_surface():
    params, rmse = calibrate_ssvi_arbitrage_free(_market_from(ARB))
    assert ssvi_is_arbitrage_free(params)


def test_plain_fit_may_be_arbitrageable():
    # The plain fit reproduces the arbitraging surface (so it is not arb-free),
    # which is exactly what the penalised variant is there to fix.
    params, _ = calibrate_ssvi(_market_from(ARB))
    assert not ssvi_is_arbitrage_free(params)


def test_arb_free_fit_on_clean_surface_matches_plain():
    market = _market_from(GOOD)
    p_plain, r_plain = calibrate_ssvi(market)
    p_af, r_af = calibrate_ssvi_arbitrage_free(market)
    assert ssvi_is_arbitrage_free(p_af)
    # A clean surface needs no penalty, so the fit quality is essentially equal.
    assert r_af == pytest.approx(r_plain, abs=1e-4)


def test_penalty_costs_some_fit_on_arb_surface():
    market = _market_from(ARB)
    _, r_plain = calibrate_ssvi(market)
    _, r_af = calibrate_ssvi_arbitrage_free(market)
    # Removing the arbitrage cannot fit the (arbitraging) market better.
    assert r_af >= r_plain - 1e-9
