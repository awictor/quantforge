"""GJR-GARCH leverage volatility model."""

import pytest

from quantforge import GJRGarchParams, gjr_garch_variance, gjr_garch_forecast


P = GJRGarchParams(omega=1e-6, alpha=0.05, beta=0.90, gamma=0.06)


def test_persistence_includes_half_gamma():
    assert P.persistence == pytest.approx(0.05 + 0.90 + 0.03)


def test_leverage_negative_shock_raises_variance_more():
    h, r = 0.0004, 0.02
    assert gjr_garch_variance(P, -r, h) > gjr_garch_variance(P, r, h)


def test_zero_gamma_is_symmetric():
    sym = GJRGarchParams(1e-6, 0.05, 0.90, 0.0)
    assert gjr_garch_variance(sym, 0.02, 0.0004) == gjr_garch_variance(sym, -0.02, 0.0004)


def test_long_run_variance_formula():
    assert P.long_run_variance == pytest.approx(1e-6 / (1 - P.persistence))


def test_forecast_reverts_to_long_run():
    lr_ann = (P.long_run_variance * 252) ** 0.5
    assert gjr_garch_forecast(P, -0.05, 0.001, 500) == pytest.approx(lr_ann, abs=1e-3)


def test_validation():
    with pytest.raises(ValueError):
        gjr_garch_forecast(P, 0.02, 0.0004, 0)
