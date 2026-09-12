"""k-nearest-neighbor classification and regression."""

import random

import pytest

from quantforge import knn_classify, knn_regress


def _two_classes(seed=1):
    rng = random.Random(seed)
    X, y = [], []
    for _ in range(50):
        X.append([rng.gauss(0, 0.5), rng.gauss(0, 0.5)])
        y.append(0)
    for _ in range(50):
        X.append([rng.gauss(5, 0.5), rng.gauss(5, 0.5)])
        y.append(1)
    return X, y


def test_k1_memorizes_training():
    X, y = _two_classes()
    assert knn_classify(X, y, X, 1) == y


def test_classifies_separable_queries():
    X, y = _two_classes()
    assert knn_classify(X, y, [[0, 0], [5, 5], [0.2, -0.1]], 5) == [0, 1, 0]


def test_regression_approximates_linear():
    X = [[float(i)] for i in range(20)]
    y = [2 * i for i in range(20)]
    preds = knn_regress(X, y, [[5.0], [10.0]], 3)
    assert abs(preds[0] - 10.0) < 1.5
    assert abs(preds[1] - 20.0) < 1.5


def test_k1_regression_is_exact_nearest():
    X = [[float(i)] for i in range(20)]
    y = [2 * i for i in range(20)]
    assert knn_regress(X, y, [[7.0]], 1) == [14.0]


def test_validation():
    X, y = _two_classes()
    with pytest.raises(ValueError):
        knn_classify(X, y, [[0, 0]], 0)
    with pytest.raises(ValueError):
        knn_classify(X, y, [[0, 0]], len(X) + 1)
    with pytest.raises(ValueError):
        knn_regress(X, [1.0], [[0, 0]], 3)      # length mismatch
