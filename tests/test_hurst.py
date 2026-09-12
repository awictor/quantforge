"""Hurst exponent via rescaled-range analysis."""

import random

import pytest

from quantforge import hurst_exponent, rescaled_range


def _white(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, 1) for _ in range(n)]


def _ar1(n, phi, seed):
    rng = random.Random(seed)
    y = [0.0]
    for _ in range(n):
        y.append(phi * y[-1] + rng.gauss(0, 1))
    return y[1:]


def test_white_noise_near_half():
    assert abs(hurst_exponent(_white(8192, 1)) - 0.5) < 0.1


def test_random_walk_near_one():
    wn = _white(8192, 1)
    walk, s = [], 0.0
    for v in wn:
        s += v
        walk.append(s)
    assert hurst_exponent(walk) > 0.8


def test_mean_reverting_below_half():
    assert hurst_exponent(_ar1(8192, -0.5, 2)) < 0.5


def test_persistent_above_half():
    assert hurst_exponent(_ar1(8192, 0.7, 3)) > 0.5


def test_rescaled_range_constant_is_zero():
    assert rescaled_range([5, 5, 5, 5]) == 0.0


def test_validation():
    with pytest.raises(ValueError):
        hurst_exponent([1.0, 2.0, 3.0])          # too short
    with pytest.raises(ValueError):
        rescaled_range([1.0])                    # window < 2
    with pytest.raises(ValueError):
        hurst_exponent(_white(100, 1), min_window=1)   # min_window < 2
