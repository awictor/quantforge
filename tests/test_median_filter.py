"""Nonlinear order-statistic filters: median, rank, Hampel."""

import random
import statistics

import pytest

from quantforge import median_filter, rank_filter, hampel_filter


def test_median_removes_single_spike():
    x = [1.0, 1.0, 1.0, 50.0, 1.0, 1.0, 1.0]
    assert median_filter(x, 3) == [1.0] * 7


def test_median_preserves_step_edge():
    step = [0.0] * 5 + [10.0] * 5
    # a median filter keeps the step sharp; a moving average would smear it
    assert median_filter(step, 3) == step


def test_median_matches_brute_force_interior():
    rng = random.Random(1)
    data = [rng.gauss(0.0, 1.0) for _ in range(50)]
    w, half = 5, 2
    mine = median_filter(data, w)
    for i in range(half, len(data) - half):
        ref = statistics.median(data[i - half:i + half + 1])
        assert abs(mine[i] - ref) < 1e-12


def test_rank_filter_min_median_max():
    r = [3.0, 1.0, 2.0, 5.0, 4.0]
    assert rank_filter(r, 3, 0) == [1.0, 1.0, 1.0, 2.0, 4.0]      # min
    assert rank_filter(r, 3, 100) == [3.0, 3.0, 5.0, 5.0, 5.0]    # max
    # interior odd windows: p50 == median
    assert rank_filter(r, 3, 50)[1:4] == median_filter(r, 3)[1:4]


def test_hampel_flags_only_the_spike():
    x = [float(i % 5) for i in range(40)]
    x[20] = 99.0
    out, idx = hampel_filter(x, window=5, n_sigmas=3.0)
    assert idx == [20]
    assert out[20] != 99.0


def test_hampel_leaves_clean_data_untouched():
    pure = [float(i % 7) for i in range(50)]
    out, idx = hampel_filter(pure, window=5, n_sigmas=3.0)
    assert out == pure
    assert idx == []


def test_validation():
    with pytest.raises(ValueError):
        median_filter([], 3)
    with pytest.raises(ValueError):
        median_filter([1.0, 2.0, 3.0], 2)        # even window
    with pytest.raises(ValueError):
        rank_filter([1.0, 2.0, 3.0], 3, 150)     # percentile out of range
    with pytest.raises(ValueError):
        hampel_filter([1.0, 2.0], n_sigmas=0.0)  # non-positive n_sigmas
