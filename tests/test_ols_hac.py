"""Robust OLS standard errors: White (HC0) and Newey-West (HAC)."""

import random

import pytest

from quantforge import white_hc0, newey_west, ols_fit


def test_coefficients_match_ols():
    rng = random.Random(3)
    X = [[rng.gauss(0, 1)] for _ in range(200)]
    y = [2 + 3 * X[t][0] + rng.gauss(0, 1) for t in range(200)]
    w = white_hc0(X, y)
    o = ols_fit(X, y)
    assert max(abs(w["coefficients"][i] - o["coefficients"][i]) for i in range(2)) < 1e-9


def test_newey_west_zero_lags_equals_white():
    rng = random.Random(3)
    X = [[rng.gauss(0, 1)] for _ in range(200)]
    y = [2 + 3 * X[t][0] + rng.gauss(0, 1) for t in range(200)]
    w = white_hc0(X, y)
    nw = newey_west(X, y, 0)
    assert max(abs(nw["std_errors"][i] - w["std_errors"][i]) for i in range(2)) < 1e-12


def test_white_reacts_to_heteroskedasticity():
    rng = random.Random(7)
    X = [[rng.uniform(0, 5)] for _ in range(300)]
    y = [1 + 2 * X[t][0] + rng.gauss(0, 1) * X[t][0] for t in range(300)]
    w = white_hc0(X, y)
    o = ols_fit(X, y)
    assert w["std_errors"][1] > o["std_errors"][1]      # robust SE larger


def test_newey_west_reacts_to_autocorrelation():
    rng = random.Random(11)
    n = 400
    X = [[rng.gauss(0, 1)] for _ in range(n)]
    e = [rng.gauss(0, 1)]
    for t in range(1, n):
        e.append(0.7 * e[t - 1] + rng.gauss(0, 1))
    y = [1 + 0.5 * X[t][0] + e[t] for t in range(n)]
    w = white_hc0(X, y)
    nw = newey_west(X, y, 8)
    assert nw["std_errors"][0] > w["std_errors"][0]     # HAC captures AR(1)


def test_cov_diagonal_positive():
    rng = random.Random(5)
    X = [[rng.gauss(0, 1), rng.gauss(0, 1)] for _ in range(150)]
    y = [1 + X[t][0] - 2 * X[t][1] + rng.gauss(0, 1) for t in range(150)]
    nw = newey_west(X, y, 4)
    assert all(nw["cov"][i][i] > 0 for i in range(3))


def test_validation():
    X = [[1.0], [2.0], [3.0]]
    y = [1.0, 2.0, 3.0]
    with pytest.raises(ValueError):
        newey_west(X, y, -1)
    with pytest.raises(ValueError):
        newey_west(X, y, 10)             # lags >= n
