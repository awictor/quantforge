"""EWMA covariance/correlation and realized beta (correlation module)."""

import random

import pytest

from quantforge import ewma_covariance, ewma_correlation, realized_beta


def _series(n=500, seed=3):
    rng = random.Random(seed)
    market = [rng.gauss(0.0, 0.01) for _ in range(n)]
    asset = [1.5 * x + rng.gauss(0.0, 0.005) for x in market]
    return asset, market


def test_realized_beta_recovers_slope():
    asset, market = _series()
    assert realized_beta(asset, market) == pytest.approx(1.5, abs=0.1)


def test_ewma_correlation_in_range():
    asset, market = _series()
    rho = ewma_correlation(asset, market)
    assert -1.0 <= rho <= 1.0
    assert rho > 0.8   # strongly positively related


def test_self_correlation_is_one():
    _, market = _series()
    assert ewma_correlation(market, market) == pytest.approx(1.0, abs=1e-9)


def test_anti_correlation_is_minus_one():
    _, market = _series()
    neg = [-x for x in market]
    assert ewma_correlation(market, neg) == pytest.approx(-1.0, abs=1e-9)


def test_ewma_covariance_self_is_positive():
    _, market = _series()
    assert ewma_covariance(market, market) > 0.0


def test_beta_of_market_to_itself_is_one():
    _, market = _series()
    assert realized_beta(market, market) == pytest.approx(1.0, abs=1e-9)


def test_validation():
    asset, market = _series(n=50)
    with pytest.raises(ValueError):
        realized_beta(asset, market[:10])
    with pytest.raises(ValueError):
        ewma_correlation(asset, market, lam=1.5)
    with pytest.raises(ValueError):
        ewma_covariance([0.01], [0.01])
