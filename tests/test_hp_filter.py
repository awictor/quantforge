"""Hodrick-Prescott trend/cycle filter."""

import random

import pytest

from quantforge import hp_filter
from quantforge.hp_filter import _solve_pentadiagonal_banded
from quantforge.portopt import _invert, _matvec


def _series(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, 1) + 0.1 * i for i in range(n)]


def _dense_solve(n, lam, y):
    D = [[0.0] * n for _ in range(n - 2)]
    for r in range(n - 2):
        D[r][r], D[r][r + 1], D[r][r + 2] = 1.0, -2.0, 1.0
    DtD = [[sum(D[k][i] * D[k][j] for k in range(n - 2)) for j in range(n)]
           for i in range(n)]
    M = [[(1.0 if i == j else 0.0) + lam * DtD[i][j] for j in range(n)]
         for i in range(n)]
    return _matvec(_invert(M), y)


def test_banded_solve_matches_dense():
    for n in (6, 12, 25):
        for lam in (1.0, 100.0, 1600.0):
            y = _series(n, n * 10 + int(lam))
            band = _solve_pentadiagonal_banded(n, lam, y)
            dense = _dense_solve(n, lam, y)
            assert max(abs(band[i] - dense[i]) for i in range(n)) < 1e-9


def test_reconstruction():
    y = _series(40, 2)
    trend, cycle = hp_filter(y, 1600)
    assert max(abs(trend[i] + cycle[i] - y[i]) for i in range(40)) < 1e-12


def test_zero_lambda_returns_data():
    y = _series(40, 3)
    trend, _ = hp_filter(y, 0.0)
    assert max(abs(trend[i] - y[i]) for i in range(40)) < 1e-9


def test_huge_lambda_gives_linear_trend():
    y = _series(40, 4)
    trend, _ = hp_filter(y, 1e12)
    # Second difference of the trend is driven to zero -> straight line.
    sd = [trend[i + 1] - 2 * trend[i] + trend[i - 1] for i in range(1, 39)]
    assert max(abs(x) for x in sd) < 1e-4


def test_higher_lambda_is_smoother():
    y = _series(60, 5)

    def curvature(t):
        return sum((t[i + 1] - 2 * t[i] + t[i - 1]) ** 2 for i in range(1, len(t) - 1))

    lo, _ = hp_filter(y, 10)
    hi, _ = hp_filter(y, 10000)
    assert curvature(hi) < curvature(lo)


def test_short_series_passthrough():
    trend, cycle = hp_filter([1.0, 2.0])
    assert trend == [1.0, 2.0]
    assert cycle == [0.0, 0.0]


def test_validation():
    with pytest.raises(ValueError):
        hp_filter([1.0, 2.0, 3.0], lam=-1.0)
