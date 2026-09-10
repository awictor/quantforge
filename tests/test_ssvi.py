"""Surface SVI (SSVI): calibration recovery and no-arbitrage conditions."""

import math

import pytest

from quantforge import (
    SSVIParams,
    ssvi_phi,
    ssvi_total_variance,
    calibrate_ssvi,
    ssvi_butterfly_free,
    ssvi_calendar_free,
    ssvi_is_arbitrage_free,
)


KS = [-0.3, -0.15, 0.0, 0.15, 0.3]
TRUE = SSVIParams(rho=-0.4, eta=1.0, gamma=0.5,
                  thetas={0.25: 0.01, 1.0: 0.045, 2.0: 0.10})


def _market(params=TRUE):
    out = []
    for t in params.thetas:
        for k in KS:
            w = params.total_variance(k, t)
            out.append((t, k, math.sqrt(w / t)))
    return out


def test_atm_total_variance_is_theta():
    # w(0, theta) = theta by construction.
    for t, theta in TRUE.thetas.items():
        assert ssvi_total_variance(0.0, theta, TRUE.rho, TRUE.eta,
                                   TRUE.gamma) == pytest.approx(theta)


def test_calibration_recovers_synthetic_surface():
    params, rmse = calibrate_ssvi(_market())
    assert rmse < 1e-4
    assert params.rho == pytest.approx(TRUE.rho, abs=1e-3)
    assert params.eta == pytest.approx(TRUE.eta, abs=1e-3)
    assert params.gamma == pytest.approx(TRUE.gamma, abs=1e-3)
    for t, theta in TRUE.thetas.items():
        assert params.thetas[t] == pytest.approx(theta, rel=1e-2)


def test_synthetic_surface_is_arbitrage_free():
    assert ssvi_is_arbitrage_free(TRUE)


def test_butterfly_condition_flags_excessive_skew():
    assert ssvi_butterfly_free(0.05, -0.4, 1.0, 0.5)
    # A huge eta blows the skew past the density-positivity bound.
    assert not ssvi_butterfly_free(0.5, -0.4, 15.0, 0.5)


def test_calendar_condition_flags_decreasing_theta():
    good = SSVIParams(rho=-0.3, eta=1.0, gamma=0.5, thetas={0.5: 0.03, 1.0: 0.06})
    bad = SSVIParams(rho=-0.3, eta=1.0, gamma=0.5, thetas={0.5: 0.06, 1.0: 0.03})
    assert ssvi_calendar_free(good)
    assert not ssvi_calendar_free(bad)


def test_phi_is_positive_and_decreasing_in_theta():
    p1 = ssvi_phi(0.02, 1.0, 0.5)
    p2 = ssvi_phi(0.10, 1.0, 0.5)
    assert p1 > 0 and p2 > 0
    assert p1 > p2  # power-law skew decays with maturity/variance


def test_negative_rho_gives_downward_skew():
    lo = TRUE.implied_vol(-0.3, 1.0)
    hi = TRUE.implied_vol(0.3, 1.0)
    assert lo > hi


def test_calibrate_requires_enough_points():
    with pytest.raises(ValueError):
        calibrate_ssvi([(1.0, 0.0, 0.2), (1.0, 0.1, 0.21)])
