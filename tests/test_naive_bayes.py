"""Gaussian naive Bayes classifier."""

import random

import pytest

from quantforge import (
    fit_gaussian_nb, predict_gaussian_nb, predict_proba_gaussian_nb,
)


def _two_classes(seed=1):
    rng = random.Random(seed)
    X, y = [], []
    for _ in range(100):
        X.append([rng.gauss(0, 1), rng.gauss(0, 1)])
        y.append("a")
    for _ in range(100):
        X.append([rng.gauss(4, 1), rng.gauss(4, 1)])
        y.append("b")
    return X, y


def test_recovers_class_means_and_priors():
    X, y = _two_classes()
    m = fit_gaussian_nb(X, y)
    assert all(abs(v) < 0.2 for v in m["means"]["a"])
    assert all(abs(v - 4) < 0.2 for v in m["means"]["b"])
    assert abs(m["priors"]["a"] - 0.5) < 1e-9


def test_classifies_and_scores_high_accuracy():
    X, y = _two_classes()
    m = fit_gaussian_nb(X, y)
    assert predict_gaussian_nb(m, [[0, 0], [4, 4]]) == ["a", "b"]
    acc = sum(1 for i in range(len(y)) if predict_gaussian_nb(m, [X[i]])[0] == y[i]) / len(y)
    assert acc > 0.95


def test_proba_sums_to_one_and_favors_nearby_class():
    X, y = _two_classes()
    m = fit_gaussian_nb(X, y)
    pr = predict_proba_gaussian_nb(m, [[0, 0]])[0]
    assert abs(sum(pr.values()) - 1.0) < 1e-12
    assert pr["a"] > 0.9


def test_validation():
    X, y = _two_classes()
    with pytest.raises(ValueError):
        fit_gaussian_nb([], [])
    with pytest.raises(ValueError):
        fit_gaussian_nb(X, ["a"])              # length mismatch
