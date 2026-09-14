"""Tests for RunningMedian, cross-checked against statistics.median after each push."""

import random
import statistics

import pytest

from quantforge.running_median import RunningMedian


def test_fuzz_continuous_stream():
    rng = random.Random(381)
    for _ in range(3000):
        n = rng.randint(1, 60)
        rm = RunningMedian()
        seen = []
        for _ in range(n):
            x = rng.uniform(-100, 100)
            rm.push(x)
            seen.append(x)
            assert abs(rm.median() - statistics.median(seen)) < 1e-9
        assert rm.count() == n
        assert len(rm) == n


def test_fuzz_integer_streams_with_dups():
    rng = random.Random(382)
    for _ in range(2000):
        n = rng.randint(1, 40)
        vals = [rng.randint(-10, 10) for _ in range(n)]
        rm = RunningMedian(vals)
        assert abs(rm.median() - statistics.median(vals)) < 1e-9


def test_even_count_averages():
    assert RunningMedian([1, 2, 3, 4]).median() == 2.5
    assert RunningMedian([5, 15, 1, 3]).median() == 4.0  # sorted 1,3,5,15 -> (3+5)/2


def test_odd_count():
    assert RunningMedian([1, 2, 3, 4, 5]).median() == 3


def test_incremental_updates():
    rm = RunningMedian()
    rm.push(10)
    assert rm.median() == 10
    rm.push(20)
    assert rm.median() == 15
    rm.push(30)
    assert rm.median() == 20


def test_constructor_from_range():
    assert RunningMedian(range(1, 6)).median() == 3


def test_single_value():
    assert RunningMedian([42]).median() == 42


def test_all_equal():
    assert RunningMedian([7, 7, 7, 7]).median() == 7


def test_descending_input():
    rm = RunningMedian()
    for x in (5, 4, 3, 2, 1):
        rm.push(x)
    assert rm.median() == 3


def test_push_returns_self():
    rm = RunningMedian()
    assert rm.push(1) is rm


def test_count_and_len():
    rm = RunningMedian([1, 2, 3])
    assert rm.count() == 3
    assert len(rm) == 3


def test_empty_median_raises():
    with pytest.raises(ValueError):
        RunningMedian().median()


def test_negative_and_mixed():
    rm = RunningMedian([-5, 0, 5, 10, -10])
    assert rm.median() == 0  # sorted -10,-5,0,5,10
