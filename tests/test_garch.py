"""Tests for GARCH(1,1) fitting and volatility forecasting."""

import math
import random

import pytest

from quantforge import fit_garch, garch_forecast, GarchParams


def _simulate_garch(omega, alpha, beta, n=4000, seed=1):
    rng = random.Random(seed)
    h = omega / (1.0 - alpha - beta)
    rets = []
    for _ in range(n):
        ret = math.sqrt(h) * rng.gauss(0, 1)
        rets.append(ret)
        h = omega + alpha * ret * ret + beta * h
    return rets


def test_recovers_persistence():
    rets = _simulate_garch(1e-6, 0.08, 0.9)
    p = fit_garch(rets)
    assert p.persistence == pytest.approx(0.98, abs=0.03)


def test_parameters_are_stationary_and_positive():
    p = fit_garch(_simulate_garch(1e-6, 0.08, 0.9))
    assert p.omega > 0
    assert p.alpha >= 0 and p.beta >= 0
    assert p.persistence < 1.0


def test_long_run_variance_matches():
    omega, alpha, beta = 1e-6, 0.08, 0.9
    p = fit_garch(_simulate_garch(omega, alpha, beta))
    true_lr = omega / (1.0 - alpha - beta)
    assert p.long_run_variance == pytest.approx(true_lr, rel=0.4)


def test_forecast_mean_reverts_to_long_run():
    p = GarchParams(omega=1e-6, alpha=0.08, beta=0.9)
    h = p.long_run_variance * 4.0   # start well above the long-run level
    lr_vol = math.sqrt(p.long_run_variance * 252)
    near = garch_forecast(p, 0.0, h, horizon=1)
    far = garch_forecast(p, 0.0, h, horizon=500)
    assert near > far                       # decays toward the long run
    assert far == pytest.approx(lr_vol, abs=1e-3)


def test_one_step_forecast_formula():
    p = GarchParams(omega=1e-6, alpha=0.08, beta=0.9)
    r_last, h_last = 0.02, 4e-4
    f = garch_forecast(p, r_last, h_last, horizon=1)
    h1 = p.omega + p.alpha * r_last ** 2 + p.beta * h_last
    assert f == pytest.approx(math.sqrt(h1 * 252))


def test_fit_rejects_short_series():
    with pytest.raises(ValueError):
        fit_garch([0.01, -0.01, 0.005])


def test_forecast_rejects_bad_horizon():
    with pytest.raises(ValueError):
        garch_forecast(GarchParams(1e-6, 0.08, 0.9), 0.0, 4e-4, horizon=0)
