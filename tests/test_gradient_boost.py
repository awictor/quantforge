"""Gradient-boosted regression trees."""

import math
import random

import pytest

from quantforge import fit_gradient_boost, predict_gradient_boost
from quantforge import fit_regression_tree, predict_regression_tree

pytestmark = pytest.mark.slow


def _mse(p, y):
    return sum((p[i] - y[i]) ** 2 for i in range(len(y))) / len(y)


def test_beats_single_tree():
    rng = random.Random(3)
    X = [[rng.uniform(0, 10)] for _ in range(300)]
    y = [math.sin(x[0]) + 0.1 * x[0] for x in X]
    gb = predict_gradient_boost(
        fit_gradient_boost(X, y, n_estimators=100, learning_rate=0.1, max_depth=3), X)
    tree = predict_regression_tree(fit_regression_tree(X, y, max_depth=3), X)
    assert _mse(gb, y) < _mse(tree, y)


def test_error_decreases_with_estimators():
    rng = random.Random(3)
    X = [[rng.uniform(0, 10)] for _ in range(300)]
    y = [math.sin(x[0]) for x in X]

    def m(n):
        return _mse(predict_gradient_boost(
            fit_gradient_boost(X, y, n_estimators=n, learning_rate=0.1, max_depth=3), X), y)

    assert m(200) < m(50) < m(10)


def test_zero_learning_rate_is_mean():
    rng = random.Random(3)
    X = [[rng.uniform(0, 10)] for _ in range(100)]
    y = [rng.gauss(0, 1) for _ in range(100)]
    p = predict_gradient_boost(fit_gradient_boost(X, y, n_estimators=50, learning_rate=0.0), X)
    assert all(abs(v - sum(y) / len(y)) < 1e-9 for v in p)


def test_constant_target():
    m = fit_gradient_boost([[1.0], [2.0], [3.0]], [5.0, 5.0, 5.0], n_estimators=10)
    assert all(abs(v - 5.0) < 1e-6 for v in predict_gradient_boost(m, [[1.5], [2.5]]))


def test_generalizes_on_linear():
    rng = random.Random(7)
    Xtr = [[rng.uniform(0, 10)] for _ in range(200)]
    ytr = [2 * x[0] + 1 for x in Xtr]
    m = fit_gradient_boost(Xtr, ytr, n_estimators=100, learning_rate=0.1, max_depth=3)
    Xte = [[rng.uniform(0, 10)] for _ in range(50)]
    yte = [2 * x[0] + 1 for x in Xte]
    assert _mse(predict_gradient_boost(m, Xte), yte) < 0.1


def test_validation():
    with pytest.raises(ValueError):
        fit_gradient_boost([[1.0]], [1.0, 2.0])
    with pytest.raises(ValueError):
        fit_gradient_boost([[1.0]], [1.0], n_estimators=0)
