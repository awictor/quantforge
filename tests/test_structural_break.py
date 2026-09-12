"""Structural-break diagnostics: CUSUM of mean and the Chow test."""

import random

import pytest

from quantforge import cusum_mean, cusum_break_detected, chow_test


def _stable(n, seed):
    rng = random.Random(seed)
    return [rng.gauss(0, 1) for _ in range(n)]


def _shift(n_each, jump, seed):
    rng = random.Random(seed)
    return ([rng.gauss(0, 1) for _ in range(n_each)]
            + [rng.gauss(jump, 1) for _ in range(n_each)])


def test_stable_series_low_false_positive_rate():
    # ~5% band; require the clear majority of stable runs to show no break.
    fp = sum(cusum_break_detected(_stable(200, s)) for s in range(40))
    assert fp <= 6


def test_level_shift_detected():
    assert cusum_break_detected(_shift(100, 5.0, 2))


def test_cusum_magnitude_larger_on_shift():
    c_shift, band = cusum_mean(_shift(100, 5.0, 2))
    c_stable, _ = cusum_mean(_stable(200, 1))
    assert max(abs(v) for v in c_shift) > band
    assert max(abs(v) for v in c_stable) < band


def test_chow_large_at_true_break():
    f, d1, d2 = chow_test(_shift(100, 5.0, 2), 100)
    assert f > 50.0
    assert d1 == 1
    assert d2 == 198


def test_chow_small_on_stable():
    f, _, _ = chow_test(_stable(200, 3), 100)
    assert f < 5.0


def test_validation():
    with pytest.raises(ValueError):
        cusum_mean([1.0, 2.0])                 # too short
    with pytest.raises(ValueError):
        cusum_mean([5.0] * 10)                 # zero variance
    with pytest.raises(ValueError):
        cusum_mean(_stable(50, 1), confidence=0.5)  # bad confidence
    with pytest.raises(ValueError):
        chow_test(_stable(50, 1), 1)           # break too close to the edge
