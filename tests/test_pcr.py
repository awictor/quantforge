"""Principal components regression."""

import random

import pytest

from quantforge import principal_components_regression, ols_fit


def test_full_components_equals_ols():
    rng = random.Random(3)
    n = 200
    X = [[rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1)] for _ in range(n)]
    y = [2 + 1.5 * X[t][0] - 0.8 * X[t][1] + 0.3 * X[t][2] + rng.gauss(0, 0.5)
         for t in range(n)]
    pcr = principal_components_regression(X, y)
    o = ols_fit(X, y)
    assert max(abs(pcr["coefficients"][i] - o["coefficients"][i + 1]) for i in range(3)) < 1e-6
    assert abs(pcr["intercept"] - o["coefficients"][0]) < 1e-6
    assert abs(pcr["explained_variance"] - 1.0) < 1e-9


def test_explained_variance_monotone():
    rng = random.Random(7)
    n = 100
    x1 = [rng.gauss(0, 1) for _ in range(n)]
    X = [[x1[t], x1[t] + rng.gauss(0, 1e-3), rng.gauss(0, 1)] for t in range(n)]
    y = [1 + 2 * x1[t] + 0.5 * X[t][2] + rng.gauss(0, 0.3) for t in range(n)]
    e1 = principal_components_regression(X, y, 1)["explained_variance"]
    e2 = principal_components_regression(X, y, 2)["explained_variance"]
    e3 = principal_components_regression(X, y, 3)["explained_variance"]
    assert e1 <= e2 <= e3
    assert abs(e3 - 1.0) < 1e-9


def test_single_predictor_matches_ols():
    rng = random.Random(5)
    X = [[rng.gauss(0, 1)] for _ in range(50)]
    y = [3 + 2 * X[t][0] + rng.gauss(0, 0.1) for t in range(50)]
    pcr = principal_components_regression(X, y)
    o = ols_fit(X, y)
    assert abs(pcr["coefficients"][0] - o["coefficients"][1]) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        principal_components_regression([[1, 2]], [1, 2])       # mismatch
    with pytest.raises(ValueError):
        principal_components_regression([[1, 2], [3, 4]], [1, 2], n_components=5)
