"""General linear-Gaussian Kalman filter and RTS smoother."""

import random

import pytest

from quantforge import kalman_filter, kalman_smoother
from quantforge.kalman import kalman_local_level


def test_reduces_to_scalar_local_level():
    rng = random.Random(1)
    obs = [rng.gauss(0, 1) for _ in range(50)]
    q, r = 0.1, 1.0
    levels, _, _ = kalman_local_level(obs, q, r)
    res = kalman_filter([[o] for o in obs], [[1.0]], [[1.0]], [[q]], [[r]],
                        [obs[0]], [[r]])
    fs = [s[0] for s in res["states"]]
    assert max(abs(fs[i] - levels[i]) for i in range(50)) < 1e-9


def test_smoother_variance_not_larger_than_filter():
    rng = random.Random(2)
    obs = [[rng.gauss(0, 1)] for _ in range(60)]
    filt = kalman_filter(obs, [[1.0]], [[1.0]], [[0.1]], [[1.0]], [0.0], [[1e6]])
    sm = kalman_smoother(obs, [[1.0]], [[1.0]], [[0.1]], [[1.0]], [0.0], [[1e6]])
    for t in range(60):
        assert sm["covariances"][t][0][0] <= filt["covariances"][t][0][0] + 1e-9


def test_constant_velocity_tracking():
    rng = random.Random(3)
    F = [[1, 1], [0, 1]]
    H = [[1, 0]]
    Q = [[1e-4, 0], [0, 1e-4]]
    R = [[0.5]]
    true = [0.5 * t for t in range(40)]
    ys = [[true[t] + rng.gauss(0, 0.5)] for t in range(40)]
    res = kalman_filter(ys, F, H, Q, R, [0.0, 0.0], [[1, 0], [0, 1]])
    est = [s[0] for s in res["states"]]
    err = sum(abs(est[t] - true[t]) for t in range(30, 40)) / 10
    assert err < 0.5
    assert abs(res["states"][-1][1] - 0.5) < 0.1     # recovered velocity


def test_log_likelihood_finite():
    rng = random.Random(4)
    obs = [[rng.gauss(0, 1)] for _ in range(30)]
    res = kalman_filter(obs, [[1.0]], [[1.0]], [[0.05]], [[1.0]], [0.0], [[1.0]])
    import math
    assert math.isfinite(res["log_likelihood"])


def test_smoother_endpoint_matches_filter():
    rng = random.Random(5)
    obs = [[rng.gauss(0, 1)] for _ in range(20)]
    filt = kalman_filter(obs, [[1.0]], [[1.0]], [[0.1]], [[1.0]], [0.0], [[1.0]])
    sm = kalman_smoother(obs, [[1.0]], [[1.0]], [[0.1]], [[1.0]], [0.0], [[1.0]])
    # The last smoothed state equals the last filtered state.
    assert abs(sm["states"][-1][0] - filt["states"][-1][0]) < 1e-12
