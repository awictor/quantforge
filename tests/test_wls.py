"""Weighted and generalized least squares."""

import random

import pytest

from quantforge import weighted_least_squares, generalized_least_squares, ols_fit


def test_equal_weights_equals_ols():
    rng = random.Random(3)
    X = [[rng.gauss(0, 1), rng.gauss(0, 1)] for _ in range(100)]
    y = [1 + 2 * X[t][0] - X[t][1] + rng.gauss(0, 1) for t in range(100)]
    w = weighted_least_squares(X, y, [1.0] * 100)
    o = ols_fit(X, y)
    assert max(abs(w["coefficients"][i] - o["coefficients"][i]) for i in range(3)) < 1e-9
    assert max(abs(w["std_errors"][i] - o["std_errors"][i]) for i in range(3)) < 1e-9


def test_gls_diagonal_equals_wls():
    rng = random.Random(9)
    n = 80
    X = [[rng.gauss(0, 1)] for _ in range(n)]
    y = [1 + 2 * X[t][0] + rng.gauss(0, 1) for t in range(n)]
    diag = [0.5 + 2 * rng.random() for _ in range(n)]
    cov = [[diag[i] if i == j else 0.0 for j in range(n)] for i in range(n)]
    g = generalized_least_squares(X, y, cov)
    w = weighted_least_squares(X, y, [1 / d for d in diag])
    assert max(abs(g["coefficients"][i] - w["coefficients"][i]) for i in range(2)) < 1e-9


def test_gls_identity_equals_ols():
    rng = random.Random(9)
    n = 60
    X = [[rng.gauss(0, 1)] for _ in range(n)]
    y = [1 + 2 * X[t][0] + rng.gauss(0, 1) for t in range(n)]
    cov = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    g = generalized_least_squares(X, y, cov)
    o = ols_fit(X, y)
    assert max(abs(g["coefficients"][i] - o["coefficients"][i]) for i in range(2)) < 1e-9


def test_gls_ar1_recovers_slope():
    rng = random.Random(11)
    n, rho = 200, 0.6
    X = [[rng.gauss(0, 1)] for _ in range(n)]
    e = [rng.gauss(0, 1)]
    for t in range(1, n):
        e.append(rho * e[t - 1] + rng.gauss(0, 1))
    y = [1 + 2 * X[t][0] + e[t] for t in range(n)]
    cov = [[rho ** abs(i - j) / (1 - rho * rho) for j in range(n)] for i in range(n)]
    g = generalized_least_squares(X, y, cov)
    assert abs(g["coefficients"][1] - 2.0) < 0.1


def test_validation():
    with pytest.raises(ValueError):
        weighted_least_squares([[1], [2]], [1, 2], [1])
    with pytest.raises(ValueError):
        weighted_least_squares([[1], [2], [3]], [1, 2, 3], [1, -1, 1])
