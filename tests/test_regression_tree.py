"""CART regression tree."""

import random

import pytest

from quantforge import fit_regression_tree, predict_regression_tree


def test_step_function_exact():
    X = [[float(i)] for i in range(10)]
    y = [0.0 if i < 5 else 10.0 for i in range(10)]
    t = fit_regression_tree(X, y, max_depth=3)
    assert predict_regression_tree(t, X) == y


def test_constant_single_leaf():
    t = fit_regression_tree([[1.0], [2.0], [3.0]], [7.0, 7.0, 7.0])
    assert t["leaf"]
    assert predict_regression_tree(t, [[1.5]]) == [7.0]


def test_deeper_lowers_error():
    rng = random.Random(3)
    X = [[rng.uniform(0, 10)] for _ in range(300)]
    y = [(x[0] - 5) ** 2 for x in X]

    def mse(depth):
        t = fit_regression_tree(X, y, max_depth=depth)
        p = predict_regression_tree(t, X)
        return sum((p[i] - y[i]) ** 2 for i in range(len(y))) / len(y)

    assert mse(8) < mse(5) < mse(1)


def test_single_split_piecewise_mean():
    X = [[0.0], [1.0], [2.0], [3.0]]
    y = [1.0, 3.0, 3.0, 5.0]
    t = fit_regression_tree(X, y, max_depth=1, min_samples=1)
    p = predict_regression_tree(t, X)
    assert p[0] == 1.0                       # left leaf = its mean
    assert abs(p[1] - 11.0 / 3) < 1e-9       # right leaf mean of 3,3,5


def test_finds_relevant_feature():
    rng = random.Random(7)
    X = [[rng.random(), rng.random()] for _ in range(200)]
    y = [10.0 if x[1] > 0.5 else 0.0 for x in X]   # depends on feature 1
    t = fit_regression_tree(X, y, max_depth=3)
    p = predict_regression_tree(t, X)
    assert sum(abs(p[i] - y[i]) for i in range(len(y))) / len(y) < 0.01


def test_validation():
    with pytest.raises(ValueError):
        fit_regression_tree([[1.0]], [1.0, 2.0])
    with pytest.raises(ValueError):
        fit_regression_tree([], [])
