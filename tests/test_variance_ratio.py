"""Lo-MacKinlay variance-ratio test."""

import random

import pytest

from quantforge import variance_ratio, variance_ratio_zstat


def _white(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, 1) for _ in range(n)]


def _ar1(n, phi, seed):
    rng = random.Random(seed)
    y = [0.0]
    for _ in range(n):
        y.append(phi * y[-1] + rng.gauss(0, 1))
    return y[1:]


def test_white_noise_ratio_near_one():
    wn = _white(5000, 1)
    assert abs(variance_ratio(wn, 2) - 1.0) < 0.1
    assert abs(variance_ratio_zstat(wn, 2)) < 3.0


def test_mean_reverting_ratio_below_one_and_rejects():
    mr = _ar1(5000, -0.4, 2)
    assert variance_ratio(mr, 2) < 1.0
    assert variance_ratio_zstat(mr, 2) < -1.96


def test_trending_ratio_above_one_and_rejects():
    tr = _ar1(5000, 0.4, 3)
    assert variance_ratio(tr, 2) > 1.0
    assert variance_ratio_zstat(tr, 2) > 1.96


def test_validation():
    wn = _white(100, 1)
    with pytest.raises(ValueError):
        variance_ratio(wn, 1)             # q < 2
    with pytest.raises(ValueError):
        variance_ratio([1.0, 2.0, 3.0], 5)  # q >= n
    with pytest.raises(ValueError):
        variance_ratio([2.0] * 50, 2)     # zero variance
