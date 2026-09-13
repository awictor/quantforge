"""Elastic-net regression (L1 + L2)."""

import random

import pytest

from quantforge import elastic_net, lasso_regression
from quantforge.ols import ols_fit


def _sparse(seed=1, n=500):
    rng = random.Random(seed)
    X = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(n)]
    y = [3.0 + 2.0 * X[i][0] - 1.5 * X[i][2] + rng.gauss(0, 0.3) for i in range(n)]
    return X, y


def test_l1_ratio_one_matches_lasso():
    X, y = _sparse()
    en = elastic_net(X, y, alpha=0.1, l1_ratio=1.0)
    la = lasso_regression(X, y, alpha=0.1)
    assert all(abs(en[k] - la[k]) < 1e-9 for k in range(len(en)))


def test_alpha_zero_matches_ols():
    X, y = _sparse()
    en = elastic_net(X, y, alpha=0.0, l1_ratio=0.5)
    ols = ols_fit(X, y)["coefficients"]
    assert all(abs(en[k] - ols[k]) < 0.05 for k in range(len(en)))


def test_still_selects_features():
    X, y = _sparse()
    en = elastic_net(X, y, alpha=0.1, l1_ratio=0.5)
    assert en[2] == 0.0 and en[4] == 0.0        # noise features zeroed
    assert abs(en[1]) > 0.5 and abs(en[3]) > 0.5


def test_grouping_effect_on_correlated_features():
    rng = random.Random(2)
    n = 500
    x1 = [rng.gauss(0, 1) for _ in range(n)]
    X = [[x1[i], x1[i] + 0.01 * rng.gauss(0, 1)] for i in range(n)]   # near-duplicates
    y = [2.0 * x1[i] + rng.gauss(0, 0.3) for i in range(n)]
    la = lasso_regression(X, y, alpha=0.1)
    en = elastic_net(X, y, alpha=0.1, l1_ratio=0.5)
    # LASSO drops one duplicate; elastic net shares weight between them.
    assert min(abs(la[1]), abs(la[2])) < 1e-6
    assert abs(en[1]) > 0.1 and abs(en[2]) > 0.1


def test_validation():
    with pytest.raises(ValueError):
        elastic_net([[1.0]], [1.0], alpha=-1.0)
    with pytest.raises(ValueError):
        elastic_net([[1.0], [2.0]], [1.0, 2.0], l1_ratio=1.5)
