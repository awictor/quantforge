"""Tests for sliding-window min/max/sum, cross-checked against a per-window brute scan."""

import random

import pytest

from quantforge.sliding_window import sliding_window_min, sliding_window_max, sliding_window_sum


def _brute(vals, k, fn):
    return [fn(vals[i:i + k]) for i in range(len(vals) - k + 1)]


def test_fuzz_vs_brute():
    rng = random.Random(391)
    for _ in range(5000):
        n = rng.randint(1, 40)
        vals = [rng.randint(-30, 30) for _ in range(n)]
        k = rng.randint(1, n)
        assert sliding_window_min(vals, k) == _brute(vals, k, min)
        assert sliding_window_max(vals, k) == _brute(vals, k, max)
        assert sliding_window_sum(vals, k) == _brute(vals, k, sum)


def test_explicit():
    v = [1, 3, -1, -3, 5, 3, 6, 7]
    assert sliding_window_max(v, 3) == [3, 3, 5, 5, 6, 7]
    assert sliding_window_min(v, 3) == [-1, -3, -3, -3, 3, 3]
    assert sliding_window_sum([1, 2, 3, 4], 2) == [3, 5, 7]


def test_window_size_one_is_identity():
    assert sliding_window_max([5, 2, 8], 1) == [5, 2, 8]
    assert sliding_window_min([5, 2, 8], 1) == [5, 2, 8]
    assert sliding_window_sum([5, 2, 8], 1) == [5, 2, 8]


def test_window_size_equals_length():
    assert sliding_window_max([4, 1, 7], 3) == [7]
    assert sliding_window_min([4, 1, 7], 3) == [1]
    assert sliding_window_sum([4, 1, 7], 3) == [12]


def test_duplicates():
    assert sliding_window_max([2, 2, 2], 2) == [2, 2]
    assert sliding_window_min([2, 2, 2], 2) == [2, 2]


def test_output_length():
    v = list(range(10))
    for k in range(1, 11):
        assert len(sliding_window_max(v, k)) == len(v) - k + 1


def test_monotone_increasing():
    v = [1, 2, 3, 4, 5]
    assert sliding_window_max(v, 2) == [2, 3, 4, 5]
    assert sliding_window_min(v, 2) == [1, 2, 3, 4]


def test_monotone_decreasing():
    v = [5, 4, 3, 2, 1]
    assert sliding_window_max(v, 2) == [5, 4, 3, 2]
    assert sliding_window_min(v, 2) == [4, 3, 2, 1]


def test_window_zero_raises():
    for fn in (sliding_window_min, sliding_window_max, sliding_window_sum):
        with pytest.raises(ValueError):
            fn([1, 2], 0)


def test_window_too_large_raises():
    for fn in (sliding_window_min, sliding_window_max, sliding_window_sum):
        with pytest.raises(ValueError):
            fn([1, 2], 5)


def test_float_values():
    v = [1.5, 2.5, 0.5, 3.5]
    assert sliding_window_max(v, 2) == [2.5, 2.5, 3.5]
    assert sliding_window_min(v, 2) == [1.5, 0.5, 0.5]
