"""Decision stump classifier."""

import random

import pytest

from quantforge import gini_impurity, fit_decision_stump, predict_decision_stump


def test_gini_bounds():
    assert gini_impurity([1, 1, 1]) == 0.0
    assert gini_impurity([0, 0, 1, 1]) == 0.5


def test_separable_on_a_single_feature():
    # Feature 1 separates the classes; feature 0 is pure noise (constant).
    X = [[0, 1], [0, 1.5], [0, 2], [0, 5], [0, 6], [0, 7]]
    y = [0, 0, 0, 1, 1, 1]
    s = fit_decision_stump(X, y)
    assert s["feature"] == 1
    assert s["gini"] == 0.0
    assert 2 < s["threshold"] < 5
    assert predict_decision_stump(s, X) == y


def test_noisy_data_reduces_gini():
    rng = random.Random(1)
    X = [[rng.gauss(0, 1)] for _ in range(100)] + [[rng.gauss(3, 1)] for _ in range(100)]
    y = [0] * 100 + [1] * 100
    s = fit_decision_stump(X, y)
    assert s["gini"] < gini_impurity(y)
    acc = sum(1 for i in range(200) if predict_decision_stump(s, [X[i]])[0] == y[i]) / 200
    assert acc > 0.9


def test_identical_rows_fall_back_to_majority():
    s = fit_decision_stump([[1], [1], [1]], [0, 0, 1])
    assert s["left_label"] == 0


def test_validation():
    with pytest.raises(ValueError):
        gini_impurity([])
    with pytest.raises(ValueError):
        fit_decision_stump([[1, 2]], [0, 1])   # length mismatch
