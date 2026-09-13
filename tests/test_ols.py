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


def test_f_equals_slope_t_squared_single_regressor():
    # With one regressor the overall F equals the slope t-statistic squared.
    X = [[1], [2], [3], [4], [5]]
    y = [2.1, 4.3, 5.9, 8.2, 9.8]
    r = ols_fit(X, y)
    assert abs(r["f_stat"] - r["t_stats"][1] ** 2) < 1e-6
    # And the F p-value equals the slope's two-sided p-value.
    assert abs(r["f_pvalue"] - r["p_values"][1]) < 1e-9


def test_pvalues_in_unit_interval_and_ordered_by_t():
    X = [[1], [2], [3], [4], [5]]
    y = [2.1, 4.3, 5.9, 8.2, 9.8]
    r = ols_fit(X, y)
    for pv in r["p_values"]:
        assert 0.0 <= pv <= 1.0
    # The larger |t| (slope) has the smaller p-value.
    assert r["p_values"][1] < r["p_values"][0]


def test_conf_int_symmetric_and_contains_coef():
    X = [[1], [2], [3], [4], [5]]
    y = [2.1, 4.3, 5.9, 8.2, 9.8]
    r = ols_fit(X, y, confidence=0.95)
    for b, (lo, hi) in zip(r["coefficients"], r["conf_int"]):
        assert lo < b < hi
        assert abs((lo + hi) / 2 - b) < 1e-9


def test_higher_confidence_widens_interval():
    X = [[1], [2], [3], [4], [5]]
    y = [2.1, 4.3, 5.9, 8.2, 9.8]
    r90 = ols_fit(X, y, confidence=0.90)
    r99 = ols_fit(X, y, confidence=0.99)
    w90 = r90["conf_int"][1][1] - r90["conf_int"][1][0]
    w99 = r99["conf_int"][1][1] - r99["conf_int"][1][0]
    assert w99 > w90


def test_significant_slope_has_tiny_pvalue():
    random.seed(1)
    X = [[i * 0.1] for i in range(60)]
    y = [3.0 + 2.0 * row[0] + random.gauss(0, 0.05) for row in X]
    r = ols_fit(X, y)
    assert r["p_values"][1] < 1e-20
    assert r["f_pvalue"] < 1e-20


def test_noise_only_slope_not_significant():
    random.seed(2)
    X = [[random.gauss(0, 1)] for _ in range(80)]
    y = [random.gauss(0, 1) for _ in range(80)]   # y independent of X
    r = ols_fit(X, y)
    assert r["p_values"][1] > 0.05


def test_validation():
    with pytest.raises(ValueError):
        ols_fit([[1.0]], [1.0])            # n <= p
    with pytest.raises(ValueError):
        ols_fit([[1.0], [2.0]], [1.0])     # X, y length mismatch
    with pytest.raises(ValueError):
        ols_fit([[1.0], [2.0], [3.0]], [1.0, 2.0, 3.0], confidence=1.5)
