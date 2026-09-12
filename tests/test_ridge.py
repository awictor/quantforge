"""Ridge (L2-penalized) regression."""

import random

import pytest

from quantforge import ridge_regression, ols_fit


def _data(seed):
    rng = random.Random(seed)
    X = [[rng.gauss(0, 1), rng.gauss(0, 1)] for _ in range(200)]
    y = [1.5 + 2.0 * r[0] - 1.0 * r[1] + rng.gauss(0, 0.5) for r in X]
    return X, y


def test_zero_alpha_matches_ols():
    X, y = _data(1)
    r = ridge_regression(X, y, 0.0)
    o = ols_fit(X, y)
    assert all(abs(r["coefficients"][i] - o["coefficients"][i]) < 1e-9 for i in range(3))


def test_larger_alpha_shrinks_slopes():
    X, y = _data(1)

    def slope_norm(alpha):
        c = ridge_regression(X, y, alpha)["coefficients"]
        return abs(c[1]) + abs(c[2])

    assert slope_norm(0.0) > slope_norm(10.0) > slope_norm(1000.0)


def test_intercept_not_penalized():
    X, y = _data(1)
    c = ridge_regression(X, y, 1e6)["coefficients"]
    assert abs(c[0] - sum(y) / len(y)) < 0.1     # intercept -> mean(y)
    assert abs(c[1]) < 0.2 and abs(c[2]) < 0.2   # slopes -> ~0


def test_handles_perfect_collinearity():
    rng = random.Random(2)
    xc = [rng.gauss(0, 1) for _ in range(100)]
    X = [[v, v] for v in xc]                     # perfectly collinear
    y = [3 * v + rng.gauss(0, 0.1) for v in xc]
    r = ridge_regression(X, y, 1.0)
    assert r["r_squared"] > 0.9
    with pytest.raises(ValueError):
        ols_fit(X, y)                            # OLS is singular here


def test_validation():
    X, y = _data(1)
    with pytest.raises(ValueError):
        ridge_regression(X, y, -1.0)
    with pytest.raises(ValueError):
        ridge_regression([[1.0]], [1.0, 2.0])    # row mismatch
