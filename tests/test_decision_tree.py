"""CART classification tree."""

import random

import pytest

from quantforge import (
    fit_decision_tree, predict_decision_tree, tree_depth,
    fit_decision_stump, predict_decision_stump,
)


def test_depth_one_matches_stump():
    X = [[1.0], [2.0], [8.0], [9.0]]
    y = [0, 0, 1, 1]
    t = fit_decision_tree(X, y, max_depth=1)
    s = fit_decision_stump(X, y)
    assert predict_decision_tree(t, X) == predict_decision_stump(s, X)


def test_solves_xor_with_depth():
    X = [[0, 0], [0, 1], [1, 0], [1, 1]]
    y = [0, 1, 1, 0]
    t = fit_decision_tree(X, y, max_depth=3, min_samples=1)
    assert predict_decision_tree(t, X) == y
    assert tree_depth(t) >= 2


def test_three_class_separable_perfect():
    rng = random.Random(1)
    X, y = [], []
    for c, (cx, cy) in enumerate([(0, 0), (10, 10), (0, 10)]):
        for _ in range(30):
            X.append([cx + rng.gauss(0, 0.5), cy + rng.gauss(0, 0.5)])
            y.append(c)
    t = fit_decision_tree(X, y, max_depth=5)
    acc = sum(1 for i in range(90) if predict_decision_tree(t, [X[i]])[0] == y[i]) / 90
    assert acc > 0.98


def test_deeper_fits_at_least_as_well():
    rng = random.Random(1)
    X, y = [], []
    for c, (cx, cy) in enumerate([(0, 0), (10, 10), (0, 10)]):
        for _ in range(30):
            X.append([cx + rng.gauss(0, 0.5), cy + rng.gauss(0, 0.5)])
            y.append(c)

    def acc(t):
        return sum(1 for i in range(90) if predict_decision_tree(t, [X[i]])[0] == y[i]) / 90

    assert acc(fit_decision_tree(X, y, max_depth=8)) >= acc(fit_decision_tree(X, y, max_depth=1))


def test_pure_node_is_leaf():
    t = fit_decision_tree([[1], [2]], [0, 0])
    assert "feature" not in t


def test_validation():
    with pytest.raises(ValueError):
        fit_decision_tree([], [])
    with pytest.raises(ValueError):
        fit_decision_tree([[1]], [0], max_depth=0)
