"""KPSS stationarity test."""

import random

import pytest

from quantforge import kpss_test


def test_white_noise_not_rejected():
    rng = random.Random(1)
    wn = [rng.gauss(0, 1) for _ in range(500)]
    _, p = kpss_test(wn)
    assert p >= 0.1                # stationary: not rejected


def test_random_walk_rejected():
    rng = random.Random(1)
    rw, acc = [], 0.0
    for _ in range(500):
        acc += rng.gauss(0, 1)
        rw.append(acc)
    eta, p = kpss_test(rw)
    assert eta > 0.463             # above the 5% level critical value
    assert p <= 0.05


def test_trend_stationary_ct_not_rejected():
    rng = random.Random(2)
    ts = [0.05 * i + rng.gauss(0, 1) for i in range(500)]
    _, p = kpss_test(ts, regression="ct")
    assert p > 0.05            # trend-stationary: not rejected at 5%


def test_trend_rejected_under_level_regression():
    # A trending series looks nonstationary to the level (constant-only) test.
    rng = random.Random(2)
    ts = [0.05 * i + rng.gauss(0, 1) for i in range(500)]
    _, p = kpss_test(ts, regression="c")
    assert p <= 0.05


def test_statistic_non_negative():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(100)]
    eta, _ = kpss_test(x)
    assert eta >= 0.0


def test_validation():
    with pytest.raises(ValueError):
        kpss_test([1.0, 2.0])
    with pytest.raises(ValueError):
        kpss_test([1.0, 2.0, 3.0, 4.0], regression="x")
