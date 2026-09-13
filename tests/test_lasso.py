"""LASSO regression by coordinate descent."""

import random

import pytest

from quantforge import lasso_regression
from quantforge.ols import ols_fit


def _sparse_data(seed=1, n=500):
    rng = random.Random(seed)
    X = [[rng.gauss(0, 1) for _ in range(4)] for _ in range(n)]
    # Only features 0 and 2 matter.
    y = [3.0 + 2.0 * X[i][0] - 1.5 * X[i][2] + rng.gauss(0, 0.3) for i in range(n)]
    return X, y


def test_alpha_zero_matches_ols():
    X, y = _sparse_data()
    b = lasso_regression(X, y, alpha=0.0)
    ols = ols_fit(X, y)["coefficients"]
    assert all(abs(b[k] - ols[k]) < 0.05 for k in range(len(b)))


def test_selects_out_noise_features():
    X, y = _sparse_data()
    b = lasso_regression(X, y, alpha=0.1)
    # Features 1 and 3 (indices 2 and 4 in b) are noise -> exactly zero.
    assert b[2] == 0.0
    assert b[4] == 0.0
    # Real features kept.
    assert abs(b[1]) > 0.5
    assert abs(b[3]) > 0.5


def test_large_alpha_zeros_all_slopes():
    X, y = _sparse_data()
    b = lasso_regression(X, y, alpha=100.0)
    assert all(abs(b[k]) < 1e-9 for k in range(1, len(b)))
    assert abs(b[0] - sum(y) / len(y)) < 1e-6      # intercept = mean(y)


def test_shrinks_monotonically():
    X, y = _sparse_data()
    b_small = lasso_regression(X, y, alpha=0.05)
    b_large = lasso_regression(X, y, alpha=0.3)
    # The real-feature slope shrinks toward zero as alpha grows.
    assert abs(b_large[1]) < abs(b_small[1])


def test_validation():
    with pytest.raises(ValueError):
        lasso_regression([[1.0]], [1.0], alpha=-1.0)
    with pytest.raises(ValueError):
        lasso_regression([], [], alpha=1.0)
    with pytest.raises(ValueError):
        lasso_regression([[1.0, 2.0], [3.0]], [1.0, 2.0])
