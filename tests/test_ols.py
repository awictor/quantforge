"""Ordinary least squares with diagnostics."""

import random

import pytest

from quantforge import ols_fit


def test_exact_line_r2_one():
    X = [[float(i)] for i in range(10)]
    y = [3 * i + 2 for i in range(10)]
    m = ols_fit(X, y)
    assert abs(m["coefficients"][0] - 2.0) < 1e-9
    assert abs(m["coefficients"][1] - 3.0) < 1e-9
    assert abs(m["r_squared"] - 1.0) < 1e-12
    assert max(abs(r) for r in m["residuals"]) < 1e-9


def test_multivariate_recovers_coefficients():
    random.seed(1)
    X = [[random.gauss(0, 1), random.gauss(0, 1)] for _ in range(200)]
    y = [1.5 + 2.0 * r[0] - 1.0 * r[1] + random.gauss(0, 0.5) for r in X]
    m = ols_fit(X, y)
    b = m["coefficients"]
    assert abs(b[0] - 1.5) < 0.15
    assert abs(b[1] - 2.0) < 0.1
    assert abs(b[2] + 1.0) < 0.1


def test_r2_bounds_and_adjustment():
    random.seed(1)
    X = [[random.gauss(0, 1), random.gauss(0, 1)] for _ in range(200)]
    y = [1.5 + 2.0 * r[0] - 1.0 * r[1] + random.gauss(0, 0.5) for r in X]
    m = ols_fit(X, y)
    assert 0.0 < m["r_squared"] < 1.0
    assert m["adj_r_squared"] <= m["r_squared"]


def test_significant_predictors_and_f():
    random.seed(1)
    X = [[random.gauss(0, 1), random.gauss(0, 1)] for _ in range(200)]
    y = [1.5 + 2.0 * r[0] - 1.0 * r[1] + random.gauss(0, 0.5) for r in X]
    m = ols_fit(X, y)
    assert all(abs(t) > 3 for t in m["t_stats"][1:])
    assert m["f_stat"] > 100


def test_simple_regression_matches_closed_form():
    random.seed(1)
    xs = [random.gauss(0, 1) for _ in range(200)]
    y = [2.0 * x + random.gauss(0, 0.5) for x in xs]
    m = ols_fit([[x] for x in xs], y)
    xm = sum(xs) / 200
    ym = sum(y) / 200
    slope = sum((xs[i] - xm) * (y[i] - ym) for i in range(200)) / sum((x - xm) ** 2 for x in xs)
    assert abs(m["coefficients"][1] - slope) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        ols_fit([[1.0]], [1.0])            # n <= p
    with pytest.raises(ValueError):
        ols_fit([[1.0], [2.0]], [1.0])     # X, y length mismatch
