"""Random forest classifier."""

import random

import pytest

from quantforge import (
    fit_random_forest, predict_random_forest,
    fit_decision_tree, predict_decision_tree,
)


def _blobs(seed=1):
    rng = random.Random(seed)
    X, y = [], []
    for c, (cx, cy) in enumerate([(0, 0), (10, 10), (0, 10)]):
        for _ in range(40):
            X.append([cx + rng.gauss(0, 0.5), cy + rng.gauss(0, 0.5)])
            y.append(c)
    return X, y


def test_separable_perfect():
    X, y = _blobs()
    f = fit_random_forest(X, y, n_trees=10, seed=42)
    acc = sum(1 for i in range(len(y)) if predict_random_forest(f, [X[i]])[0] == y[i]) / len(y)
    assert acc > 0.98


def test_reproducible():
    X, y = _blobs()
    fa = fit_random_forest(X, y, n_trees=10, seed=7)
    fb = fit_random_forest(X, y, n_trees=10, seed=7)
    assert predict_random_forest(fa, X) == predict_random_forest(fb, X)


def test_generalizes_at_least_as_well_as_single_tree():
    rng = random.Random(2)

    def gen(n):
        X, y = [], []
        for _ in range(n):
            lab = rng.randint(0, 1)
            cx = 2 if lab else 0
            X.append([cx + rng.gauss(0, 1.5), rng.gauss(0, 1.5)])
            y.append(lab)
        return X, y

    Xtr, ytr = gen(200)
    Xte, yte = gen(200)
    forest = fit_random_forest(Xtr, ytr, n_trees=25, max_depth=6, seed=1)
    tree = fit_decision_tree(Xtr, ytr, max_depth=6)
    af = sum(1 for i in range(200) if predict_random_forest(forest, [Xte[i]])[0] == yte[i]) / 200
    at = sum(1 for i in range(200) if predict_decision_tree(tree, [Xte[i]])[0] == yte[i]) / 200
    assert af >= at - 0.02


def test_single_tree_forest():
    X, y = _blobs()
    assert len(fit_random_forest(X, y, n_trees=1)["trees"]) == 1


def test_validation():
    X, y = _blobs()
    with pytest.raises(ValueError):
        fit_random_forest(X, y, n_trees=0)
    with pytest.raises(ValueError):
        fit_random_forest([], [])
