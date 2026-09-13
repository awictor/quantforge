"""Treynor, Jensen's alpha, M-squared, market beta."""

import math
import random

import pytest

from quantforge import (market_beta, treynor_ratio, jensens_alpha, m_squared,
                        sharpe_ratio, annualized_volatility)


def _data(seed=1, n=1000):
    rng = random.Random(seed)
    mkt = [rng.gauss(0.0004, 0.01) for _ in range(n)]
    port = [1.5 * mkt[i] + rng.gauss(0, 0.002) for i in range(n)]
    return port, mkt


def test_self_beta_is_one():
    _, mkt = _data()
    assert abs(market_beta(mkt, mkt) - 1.0) < 1e-9


def test_beta_recovers_leverage():
    port, mkt = _data()
    assert abs(market_beta(port, mkt) - 1.5) < 0.1


def test_treynor_matches_manual():
    port, mkt = _data()
    b = market_beta(port, mkt)
    ann = (sum(port) / len(port)) * 252
    assert abs(treynor_ratio(port, mkt) - ann / b) < 1e-9


def test_jensen_alpha_near_zero_for_pure_beta():
    port, mkt = _data()
    assert abs(jensens_alpha(port, mkt)) < 0.05


def test_jensen_alpha_detects_added_alpha():
    _, mkt = _data()
    port = [1.0 * mkt[i] + 0.001 for i in range(len(mkt))]   # +0.1%/period
    assert abs(jensens_alpha(port, mkt) - 0.001 * 252) < 0.02


def test_m_squared_equals_sharpe_times_market_vol():
    port, mkt = _data()
    mv = annualized_volatility(mkt, 252)
    assert abs(m_squared(port, mkt) - sharpe_ratio(port) * mv) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        market_beta([1.0], [1.0])
    with pytest.raises(ValueError):
        treynor_ratio([-0.01, 0.01], [0.01, -0.01])    # negative beta
    with pytest.raises(ValueError):
        m_squared([0.01, 0.01, 0.01], [0.01, 0.02, 0.03])   # zero-variance port
