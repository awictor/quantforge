"""EGARCH log-variance volatility model."""

import math

import pytest

from quantforge import EGarchParams, egarch_variance, egarch_forecast


P = EGarchParams(omega=-0.1, alpha=0.15, beta=0.95, gamma=-0.08)


def test_variance_always_positive():
    assert egarch_variance(P, 0.02, 0.0004) > 0
    # Even absurd parameters yield a positive variance (log formulation).
    assert egarch_variance(EGarchParams(-5, 2, -0.5, 3), 0.05, 0.0004) > 0


def test_leverage_negative_shock_raises_variance():
    h = 0.0004
    assert egarch_variance(P, -0.02, h) > egarch_variance(P, 0.02, h)


def test_zero_gamma_symmetric():
    sym = EGarchParams(-0.1, 0.15, 0.95, 0.0)
    assert egarch_variance(sym, 0.02, 0.0004) == egarch_variance(sym, -0.02, 0.0004)


def test_forecast_reverts_to_unconditional():
    log_lr = P.omega / (1 - P.beta)
    lr_ann = math.sqrt(math.exp(log_lr) * 252)
    assert egarch_forecast(P, -0.05, 0.001, 500) == pytest.approx(lr_ann, abs=1e-3)


def test_validation():
    with pytest.raises(ValueError):
        egarch_variance(P, 0.02, 0)
    with pytest.raises(ValueError):
        egarch_forecast(EGarchParams(-0.1, 0.15, 1.0, 0), 0.02, 0.0004)
