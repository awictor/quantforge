"""Poisson regression (GLM, log link, IRLS)."""

import math
import random

import pytest

from quantforge import poisson_regression, poisson_predict


def _sample_poisson(lam, rng):
    k = 0
    p = math.exp(-lam)
    s = p
    u = rng.random()
    while u > s and k < 100000:
        k += 1
        p *= lam / k
        s += p
    return k


def test_recovers_log_linear_rate():
    rng = random.Random(3)
    X = [[rng.uniform(-1, 1)] for _ in range(2000)]
    y = [_sample_poisson(math.exp(0.5 + 0.8 * x[0]), rng) for x in X]
    m = poisson_regression(X, y)
    assert abs(m["coefficients"][0] - 0.5) < 0.1
    assert abs(m["coefficients"][1] - 0.8) < 0.1


def test_intercept_only_is_log_mean():
    y = [3, 5, 2, 4, 6, 3, 4, 5]
    m = poisson_regression([[] for _ in y], y)
    assert abs(math.exp(m["coefficients"][0]) - sum(y) / len(y)) < 1e-6


def test_predictions_positive():
    rng = random.Random(5)
    X = [[rng.uniform(-2, 2)] for _ in range(50)]
    y = [_sample_poisson(math.exp(0.3 + 0.5 * x[0]), rng) for x in X]
    m = poisson_regression(X, y)
    assert all(p > 0 for p in poisson_predict(m, [[-1.0], [0.0], [1.0]]))


def test_full_beats_null_loglik():
    rng = random.Random(7)
    X = [[rng.uniform(-2, 2)] for _ in range(80)]
    y = [_sample_poisson(math.exp(0.3 + 0.6 * x[0]), rng) for x in X]
    full = poisson_regression(X, y)
    null = poisson_regression([[] for _ in y], y)
    assert full["log_likelihood"] >= null["log_likelihood"] - 1e-6


def test_validation():
    with pytest.raises(ValueError):
        poisson_regression([[1.0]], [1, 2])
    with pytest.raises(ValueError):
        poisson_regression([[1.0]], [-1])
