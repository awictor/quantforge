"""Tracking error and information ratio (perfmetrics module)."""

import math
import random
import statistics

import pytest

from quantforge import tracking_error, information_ratio, sharpe_ratio


def _series(n=300, seed=3):
    rng = random.Random(seed)
    bench = [rng.gauss(0.0, 0.01) for _ in range(n)]
    port = [b + rng.gauss(0.0002, 0.003) for b in bench]
    return port, bench


def test_tracking_error_matches_active_stdev():
    r, b = _series()
    active = [x - y for x, y in zip(r, b)]
    assert tracking_error(r, b) == pytest.approx(
        statistics.stdev(active) * math.sqrt(252), abs=1e-9)


def test_information_ratio_is_active_sharpe():
    r, b = _series()
    active = [x - y for x, y in zip(r, b)]
    assert information_ratio(r, b) == pytest.approx(sharpe_ratio(active), abs=1e-9)


def test_positive_active_gives_positive_ir():
    r, b = _series()
    assert information_ratio(r, b) > 0.0


def test_zero_active_tracking_error_zero():
    _, b = _series()
    assert tracking_error(b, b) == 0.0


def test_validation():
    r, b = _series(n=50)
    with pytest.raises(ValueError):
        information_ratio(b, b)          # zero active variance
    with pytest.raises(ValueError):
        tracking_error(r, b[:10])        # length mismatch
    with pytest.raises(ValueError):
        information_ratio([0.01], [0.0])  # need >= 2
