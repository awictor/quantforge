"""Linear quantile regression by IRLS."""

import random

import pytest

from quantforge import quantile_regression
from quantforge.prob_forecast import pinball_loss
from quantforge.ols import ols_fit


def _data(seed=1, n=2000):
    rng = random.Random(seed)
    X = [[rng.gauss(0, 1)] for _ in range(n)]
    y = [2.0 + 3.0 * X[i][0] + rng.gauss(0, 1) for i in range(n)]
    return X, y


def test_median_recovers_slope_and_intercept():
    X, y = _data()
    b = quantile_regression(X, y, 0.5)
    assert abs(b[0] - 2.0) < 0.15
    assert abs(b[1] - 3.0) < 0.1


def test_intercept_increases_with_tau():
    X, y = _data()
    b10 = quantile_regression(X, y, 0.1)
    b50 = quantile_regression(X, y, 0.5)
    b90 = quantile_regression(X, y, 0.9)
    assert b10[0] < b50[0] < b90[0]
    # Slope stays near 3 across quantiles (homoskedastic noise).
    for b in (b10, b50, b90):
        assert abs(b[1] - 3.0) < 0.15


def test_residual_quantile_property():
    X, y = _data()
    n = len(y)
    for tau in (0.1, 0.5, 0.9):
        b = quantile_regression(X, y, tau)
        frac = sum(1 for i in range(n) if y[i] < b[0] + b[1] * X[i][0]) / n
        assert abs(frac - tau) < 0.03


def test_beats_ols_on_pinball_loss():
    X, y = _data()
    n = len(y)
    tau = 0.9
    bq = quantile_regression(X, y, tau)
    bo = ols_fit(X, y)["coefficients"]
    fq = [bq[0] + bq[1] * X[i][0] for i in range(n)]
    fo = [bo[0] + bo[1] * X[i][0] for i in range(n)]
    assert pinball_loss(y, fq, tau) <= pinball_loss(y, fo, tau)


def test_multivariate():
    rng = random.Random(3)
    n = 1500
    X = [[rng.gauss(0, 1), rng.gauss(0, 1)] for _ in range(n)]
    y = [1.0 + 2.0 * X[i][0] - 1.0 * X[i][1] + rng.gauss(0, 0.5) for i in range(n)]
    b = quantile_regression(X, y, 0.5)
    assert abs(b[1] - 2.0) < 0.1
    assert abs(b[2] + 1.0) < 0.1


def test_validation():
    with pytest.raises(ValueError):
        quantile_regression([[1.0]], [1.0], tau=0.5)      # n <= p
    with pytest.raises(ValueError):
        quantile_regression([[1.0], [2.0], [3.0]], [1.0, 2.0, 3.0], tau=1.0)
