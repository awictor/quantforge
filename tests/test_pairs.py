"""Pairs trading: hedge ratio, spread, OU half-life."""

import random

import pytest

from quantforge import pairs_hedge_ratio, spread_series, ou_half_life, spread_zscore


def _cointegrated():
    random.seed(6)
    x, v = [], 100.0
    for _ in range(500):
        v += random.gauss(0, 1)
        x.append(v)
    s, noise = 0.0, []
    for _ in range(500):
        s = 0.9 * s + random.gauss(0, 0.5)   # AR(1), kappa ~ 0.1
        noise.append(s)
    y = [2.0 * x[i] + noise[i] for i in range(500)]
    return y, x


def test_pairs_hedge_ratio_recovers_beta():
    y, x = _cointegrated()
    assert pairs_hedge_ratio(y, x) == pytest.approx(2.0, abs=0.05)


def test_half_life_in_range():
    y, x = _cointegrated()
    sp = spread_series(y, x)
    hl = ou_half_life(sp)
    assert 3 < hl < 15   # ln(2)/0.1 ~ 6.9


def test_spread_zscore_finite():
    y, x = _cointegrated()
    z = spread_zscore(spread_series(y, x))
    assert isinstance(z, float)


def test_constant_spread_not_mean_reverting():
    with pytest.raises(ValueError):
        ou_half_life([5.0, 5.0, 5.0, 5.0])


def test_validation():
    with pytest.raises(ValueError):
        pairs_hedge_ratio([1, 2], [1])
    with pytest.raises(ValueError):
        spread_zscore([1.0])
