"""Time-varying regression slope via Kalman filter."""

import random

import pytest

from quantforge import kalman_regression_beta


def test_static_slope_equals_ols():
    # Q = 0 with a diffuse prior is recursive least squares -> OLS slope.
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(200)]
    y = [2.3 * xi + rng.gauss(0, 0.5) for xi in x]
    betas, _ = kalman_regression_beta(x, y, 0.0, 0.5, beta0=0.0, p0=1e9)
    ols = sum(x[i] * y[i] for i in range(200)) / sum(xi * xi for xi in x)
    assert abs(betas[-1] - ols) < 1e-6


def test_posterior_variance_shrinks():
    rng = random.Random(4)
    x = [rng.gauss(0, 1) for _ in range(200)]
    y = [1.5 * xi + rng.gauss(0, 0.5) for xi in x]
    _, var = kalman_regression_beta(x, y, 0.0, 0.5, p0=1e9)
    assert var[-1] < var[0]


def test_tracks_regime_switch():
    rng = random.Random(2)
    xs = [rng.gauss(0, 1) for _ in range(400)]
    ys = [(1.0 if i < 200 else 3.0) * xi + rng.gauss(0, 0.3)
          for i, xi in enumerate(xs)]
    betas, _ = kalman_regression_beta(xs, ys, 0.02, 0.09, beta0=1.0, p0=1.0)
    early = sum(betas[150:200]) / 50
    late = sum(betas[350:400]) / 50
    assert abs(early - 1.0) < 0.25
    assert abs(late - 3.0) < 0.4


def test_static_filter_does_not_track_regime():
    rng = random.Random(2)
    xs = [rng.gauss(0, 1) for _ in range(400)]
    ys = [(1.0 if i < 200 else 3.0) * xi + rng.gauss(0, 0.3)
          for i, xi in enumerate(xs)]
    betas, _ = kalman_regression_beta(xs, ys, 0.0, 0.09, beta0=1.0, p0=1e9)
    assert 1.0 < betas[-1] < 3.0   # settles at a blend, not the late regime


def test_validation():
    with pytest.raises(ValueError):
        kalman_regression_beta([1.0, 2.0], [1.0], 0.0, 1.0)
    with pytest.raises(ValueError):
        kalman_regression_beta([1.0], [1.0], -1.0, 1.0)
    with pytest.raises(ValueError):
        kalman_regression_beta([1.0], [1.0], 0.0, 0.0)
