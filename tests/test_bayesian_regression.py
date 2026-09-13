"""Bayesian linear regression (conjugate Gaussian posterior)."""

import random

import pytest

from quantforge import (bayesian_linear_regression, bayesian_predict,
                        ridge_regression, ols_fit)


def _data(seed=3, n=100):
    rng = random.Random(seed)
    X = [[rng.gauss(0, 1), rng.gauss(0, 1)] for _ in range(n)]
    y = [2 + 3 * X[t][0] - X[t][1] + rng.gauss(0, 0.5) for t in range(n)]
    return X, y


def test_weak_prior_approaches_ols():
    X, y = _data()
    m = bayesian_linear_regression(X, y, alpha=1e-6, beta_noise=1.0)
    o = ols_fit(X, y)["coefficients"]
    assert max(abs(m["mean"][i] - o[i]) for i in range(3)) < 1e-5


def test_posterior_std_shrinks_with_data():
    X, y = _data(n=100)
    m50 = bayesian_linear_regression(X[:50], y[:50], alpha=1.0, beta_noise=4.0)
    m100 = bayesian_linear_regression(X, y, alpha=1.0, beta_noise=4.0)
    assert m100["std"][1] < m50["std"][1]


def test_strong_prior_shrinks_coefficients():
    X, y = _data()
    mw = bayesian_linear_regression(X, y, alpha=0.01, beta_noise=1.0)
    ms = bayesian_linear_regression(X, y, alpha=100.0, beta_noise=1.0)
    assert abs(ms["mean"][1]) < abs(mw["mean"][1])


def test_predictive_variance_widens_on_extrapolation():
    m = bayesian_linear_regression([[float(i)] for i in range(20)],
                                   [2 * i + 1 for i in range(20)],
                                   alpha=0.01, beta_noise=1.0)
    _, var_near = bayesian_predict(m, [10.0])
    _, var_far = bayesian_predict(m, [100.0])
    assert var_far > var_near
    assert var_near >= 1.0                    # at least the noise variance


def test_matches_ridge_with_penalized_intercept():
    # With the intercept in X (no auto-intercept), the posterior mean matches ridge.
    rng = random.Random(5)
    X = [[1.0, rng.gauss(0, 1)] for _ in range(80)]     # explicit intercept column
    y = [2 + 3 * X[t][1] + rng.gauss(0, 0.5) for t in range(80)]
    m = bayesian_linear_regression(X, y, alpha=2.0, beta_noise=4.0, add_intercept=False)
    r = ridge_regression(X, y, alpha=2.0 / 4.0, add_intercept=False)["coefficients"]
    assert max(abs(m["mean"][i] - r[i]) for i in range(2)) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        bayesian_linear_regression([[1.0]], [1.0, 2.0])
    with pytest.raises(ValueError):
        bayesian_linear_regression([[1.0]], [1.0], alpha=0)
