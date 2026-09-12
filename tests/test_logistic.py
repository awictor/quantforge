"""Logistic regression by IRLS."""

import random

import pytest

from quantforge import fit_logistic, predict_proba
from quantforge.logistic import _sigmoid


def _dataset(seed, b=(0.5, 2.0, -1.0), n=3000):
    rng = random.Random(seed)
    X, y = [], []
    for _ in range(n):
        x1, x2 = rng.gauss(0, 1), rng.gauss(0, 1)
        p = _sigmoid(b[0] + b[1] * x1 + b[2] * x2)
        y.append(1 if rng.random() < p else 0)
        X.append([x1, x2])
    return X, y


def test_recovers_coefficients():
    X, y = _dataset(1)
    m = fit_logistic(X, y)
    b = m["coefficients"]
    assert m["converged"]
    assert abs(b[0] - 0.5) < 0.2
    assert abs(b[1] - 2.0) < 0.25
    assert abs(b[2] + 1.0) < 0.2


def test_probabilities_in_unit_interval():
    X, y = _dataset(1)
    m = fit_logistic(X, y)
    assert all(0.0 < p < 1.0 for p in predict_proba(m, X))


def test_accuracy_above_baseline():
    X, y = _dataset(1)
    m = fit_logistic(X, y)
    preds = [1 if p > 0.5 else 0 for p in predict_proba(m, X)]
    acc = sum(1 for i in range(len(y)) if preds[i] == y[i]) / len(y)
    assert acc > 0.75


def test_separable_data_classified_perfectly():
    rng = random.Random(2)
    X = [[rng.gauss(0, 1)] for _ in range(200)]
    y = [1 if row[0] > 0 else 0 for row in X]
    m = fit_logistic(X, y)
    p = predict_proba(m, X)
    acc = sum(1 for i in range(200) if (p[i] > 0.5) == (y[i] == 1)) / 200
    assert acc > 0.98


def test_positive_coefficient_is_monotone():
    X, y = _dataset(1)
    m = fit_logistic(X, y)
    assert predict_proba(m, [[2.0, 0.0]])[0] > predict_proba(m, [[-2.0, 0.0]])[0]


def test_validation():
    X, y = _dataset(1, n=50)
    with pytest.raises(ValueError):
        fit_logistic(X, [0.5] * len(X))       # non-binary
    with pytest.raises(ValueError):
        fit_logistic([[1.0]], [1, 0])         # row mismatch
