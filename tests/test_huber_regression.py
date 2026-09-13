"""Huber robust regression (IRLS M-estimator)."""

import random

import pytest

from quantforge import huber_regression, ols_fit


def test_clean_data_near_ols():
    rng = random.Random(3)
    X = [[rng.gauss(0, 1)] for _ in range(200)]
    y = [2 + 3 * X[t][0] + rng.gauss(0, 0.5) for t in range(200)]
    h = huber_regression(X, y)
    o = ols_fit(X, y)
    assert max(abs(h["coefficients"][i] - o["coefficients"][i]) for i in range(2)) < 0.1


def test_resists_outliers():
    rng = random.Random(7)
    X = [[float(i)] for i in range(50)]
    y = [2 * X[i][0] + 1 + rng.gauss(0, 0.3) for i in range(50)]
    for i in [5, 15, 25, 35]:
        y[i] += 200                       # vertical outliers
    h = huber_regression(X, y)
    o = ols_fit(X, y)
    assert abs(h["coefficients"][1] - 2) < abs(o["coefficients"][1] - 2)
    assert abs(h["coefficients"][1] - 2) < 0.1


def test_large_delta_approaches_ols():
    rng = random.Random(5)
    X = [[rng.gauss(0, 1)] for _ in range(100)]
    y = [1 + 2 * X[t][0] + rng.gauss(0, 1) for t in range(100)]
    h = huber_regression(X, y, delta=100.0)
    o = ols_fit(X, y)
    assert max(abs(h["coefficients"][i] - o["coefficients"][i]) for i in range(2)) < 1e-3


def test_exact_line():
    X = [[float(i)] for i in range(10)]
    y = [3 * X[i][0] - 2 for i in range(10)]
    h = huber_regression(X, y)
    assert abs(h["coefficients"][0] + 2) < 1e-6 and abs(h["coefficients"][1] - 3) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        huber_regression([[1.0]], [1.0, 2.0])
    with pytest.raises(ValueError):
        huber_regression([[1.0]], [1.0], delta=0)
