"""Tests for binary/ternary search primitives, cross-checked against brute scans."""

import random

import pytest

from quantforge.binary_search import (
    first_true,
    last_true,
    ternary_search_int_max,
    ternary_search_int_min,
)


def test_fuzz_first_true():
    rng = random.Random(201)
    for _ in range(5000):
        lo = rng.randint(-20, 20)
        hi = lo + rng.randint(0, 40)
        thr = rng.randint(lo - 2, hi + 2)
        pred = lambda x: x >= thr
        brute = next((x for x in range(lo, hi + 1) if pred(x)), hi + 1)
        assert first_true(lo, hi, pred) == brute


def test_fuzz_last_true():
    rng = random.Random(202)
    for _ in range(5000):
        lo = rng.randint(-20, 20)
        hi = lo + rng.randint(0, 40)
        thr = rng.randint(lo - 2, hi + 2)
        pred = lambda x: x <= thr
        brute = lo - 1
        for x in range(lo, hi + 1):
            if pred(x):
                brute = x
        assert last_true(lo, hi, pred) == brute


def test_fuzz_ternary_max():
    rng = random.Random(203)
    for _ in range(5000):
        lo = rng.randint(-30, 30)
        hi = lo + rng.randint(0, 60)
        peak = rng.randint(lo, hi)
        f = lambda x: -((x - peak) ** 2)
        got = ternary_search_int_max(lo, hi, f)
        brute = max(range(lo, hi + 1), key=f)
        assert f(got) == f(brute)


def test_fuzz_ternary_min():
    rng = random.Random(204)
    for _ in range(5000):
        lo = rng.randint(-30, 30)
        hi = lo + rng.randint(0, 60)
        valley = rng.randint(lo, hi)
        f = lambda x: (x - valley) ** 2
        got = ternary_search_int_min(lo, hi, f)
        brute = min(range(lo, hi + 1), key=f)
        assert f(got) == f(brute)


def test_binary_search_the_answer():
    assert first_true(0, 100, lambda x: x * x >= 50) == 8
    assert last_true(0, 100, lambda x: x * x <= 50) == 7


def test_all_false_and_all_true():
    assert first_true(0, 10, lambda x: False) == 11
    assert first_true(0, 10, lambda x: True) == 0
    assert last_true(0, 10, lambda x: False) == -1
    assert last_true(0, 10, lambda x: True) == 10


def test_ternary_peak_at_boundary():
    assert ternary_search_int_max(0, 10, lambda x: x) == 10
    assert ternary_search_int_max(0, 10, lambda x: -x) == 0
    assert ternary_search_int_min(0, 10, lambda x: x) == 0
    assert ternary_search_int_min(0, 10, lambda x: -x) == 10


def test_single_point():
    assert first_true(5, 5, lambda x: True) == 5
    assert first_true(5, 5, lambda x: False) == 6
    assert last_true(5, 5, lambda x: True) == 5
    assert ternary_search_int_max(5, 5, lambda x: x) == 5
    assert ternary_search_int_min(5, 5, lambda x: x) == 5


def test_negative_range():
    assert first_true(-100, -1, lambda x: x >= -10) == -10
    assert last_true(-100, -1, lambda x: x <= -10) == -10


def test_lo_greater_than_hi_raises():
    for fn in (first_true, last_true, ternary_search_int_max, ternary_search_int_min):
        with pytest.raises(ValueError):
            fn(5, 3, lambda x: True)


def test_first_true_call_count_logarithmic():
    calls = []

    def pred(x):
        calls.append(x)
        return x >= 1_000_000

    first_true(0, 1_000_000_000, pred)
    # log2(1e9) ~ 30; well under a linear scan
    assert len(calls) < 40
