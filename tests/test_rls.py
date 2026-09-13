"""Recursive least squares."""

import random

import pytest

from quantforge import RecursiveLeastSquares, recursive_least_squares, ols_fit


def test_matches_batch_ols():
    rng = random.Random(3)
    n = 300
    X = [[1.0, rng.gauss(0, 1), rng.gauss(0, 1)] for _ in range(n)]
    y = [2 + 3 * X[t][1] - 1.5 * X[t][2] + rng.gauss(0, 0.5) for t in range(n)]
    beta = recursive_least_squares(X, y, forgetting=1.0, delta=1e7)
    o = ols_fit([[row[1], row[2]] for row in X], y, add_intercept=True)
    assert max(abs(beta[i] - o["coefficients"][i]) for i in range(3)) < 1e-6


def test_exact_fit_through_points():
    X = [[1.0, 1.0], [1.0, 2.0], [1.0, 3.0]]
    y = [3, 5, 7]                            # y = 1 + 2x
    beta = recursive_least_squares(X, y, delta=1e9)
    assert abs(beta[0] - 1.0) < 1e-6 and abs(beta[1] - 2.0) < 1e-6


def test_predict():
    rng = random.Random(5)
    rls = RecursiveLeastSquares(2, delta=1e7)
    for _ in range(200):
        x = [1.0, rng.gauss(0, 1)]
        rls.update(x, 5 + 2 * x[1] + rng.gauss(0, 0.1))
    assert abs(rls.predict([1.0, 1.0]) - 7.0) < 0.1


def test_forgetting_tracks_regime_change():
    rng = random.Random(7)
    rls = RecursiveLeastSquares(2, forgetting=0.95, delta=1e3)
    for _ in range(150):
        x = [1.0, rng.gauss(0, 1)]
        rls.update(x, 1 + 2 * x[1] + rng.gauss(0, 0.1))
    for _ in range(150):
        x = [1.0, rng.gauss(0, 1)]
        rls.update(x, 1 - 3 * x[1] + rng.gauss(0, 0.1))
    assert abs(rls.beta[1] - (-3.0)) < 0.3     # adapts to the new slope


def test_validation():
    with pytest.raises(ValueError):
        RecursiveLeastSquares(0)
    with pytest.raises(ValueError):
        RecursiveLeastSquares(2, forgetting=1.5)
    with pytest.raises(ValueError):
        RecursiveLeastSquares(2).update([1.0], 3.0)
