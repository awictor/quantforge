"""Scalar local-level Kalman filter."""

import random

import pytest

from quantforge import kalman_local_level, kalman_steady_state_gain


def test_gain_converges_to_steady_state():
    Q, R = 0.01, 1.0
    _, _, gains = kalman_local_level([0.0] * 5000, Q, R)
    assert abs(gains[-1] - kalman_steady_state_gain(Q, R)) < 1e-6


def test_perfect_observations_track_data():
    # R -> 0: gain -> 1, the filter reproduces each observation.
    random.seed(1)
    y = [random.gauss(0, 1) for _ in range(50)]
    levels, _, gains = kalman_local_level(y, 1.0, 1e-9)
    assert max(abs(levels[i] - y[i]) for i in range(1, 50)) < 1e-4
    assert gains[-1] > 0.999


def test_zero_process_noise_is_running_mean():
    # Q = 0 with a diffuse prior: the estimate is the running mean and the gain
    # decays as 1/(t+1) (recursive least squares).
    random.seed(1)
    y = [random.gauss(0, 1) for _ in range(50)]
    levels, _, gains = kalman_local_level(y, 0.0, 1.0, x0=0.0, p0=1e12)
    run_mean = [sum(y[: i + 1]) / (i + 1) for i in range(50)]
    assert abs(gains[5] - 1.0 / 6) < 1e-6
    assert max(abs(levels[i] - run_mean[i]) for i in range(50)) < 1e-4


def test_variance_positive_and_decreasing_early():
    _, variances, _ = kalman_local_level([0.0] * 50, 0.01, 1.0)
    assert all(v > 0 for v in variances)
    assert variances[10] < variances[1]


def test_gain_in_unit_interval():
    random.seed(3)
    y = [random.gauss(0, 1) for _ in range(100)]
    _, _, gains = kalman_local_level(y, 0.05, 1.0)
    assert all(0.0 < g < 1.0 for g in gains)


def test_steady_state_gain_monotone_in_snr():
    # Higher process-to-observation noise ratio -> higher gain.
    g_lo = kalman_steady_state_gain(0.01, 1.0)
    g_hi = kalman_steady_state_gain(1.0, 1.0)
    assert 0.0 < g_lo < g_hi < 1.0


def test_empty_and_validation():
    assert kalman_local_level([], 0.1, 1.0) == ([], [], [])
    with pytest.raises(ValueError):
        kalman_local_level([1.0], 0.1, 0.0)      # R = 0
    with pytest.raises(ValueError):
        kalman_local_level([1.0], -1.0, 1.0)     # Q < 0
    with pytest.raises(ValueError):
        kalman_steady_state_gain(0.1, 0.0)
