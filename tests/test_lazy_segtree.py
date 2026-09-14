"""Tests for the lazy-propagation segment tree, cross-checked against a brute array."""

import random

import pytest

from quantforge.lazy_segtree import LazySegmentTree


def _ident(mode):
    return 0 if mode == "sum" else (float("inf") if mode == "min" else float("-inf"))


def _agg(mode, xs):
    if not xs:
        return _ident(mode)
    return {"sum": sum, "min": min, "max": max}[mode](xs)


def test_fuzz_vs_brute():
    rng = random.Random(311)
    for mode in ("sum", "min", "max"):
        for _ in range(3000):
            n = rng.randint(1, 40)
            arr = [rng.randint(-20, 20) for _ in range(n)]
            st = LazySegmentTree(arr, mode)
            for _ in range(20):
                if rng.random() < 0.5:
                    lo = rng.randint(0, n)
                    hi = rng.randint(lo, n)
                    d = rng.randint(-10, 10)
                    st.update(lo, hi, d)
                    for i in range(lo, hi):
                        arr[i] += d
                else:
                    lo = rng.randint(0, n)
                    hi = rng.randint(lo, n)
                    assert st.query(lo, hi) == _agg(mode, arr[lo:hi])


def test_sum_explicit():
    st = LazySegmentTree([1, 2, 3, 4, 5], "sum")
    assert st.query(0, 5) == 15
    assert st.query(1, 4) == 9
    st.update(1, 4, 10)  # [1, 12, 13, 14, 5]
    assert st.query(0, 5) == 45
    assert st.query(1, 4) == 39


def test_min_explicit():
    st = LazySegmentTree([5, 3, 8, 1, 9], "min")
    assert st.query(0, 5) == 1
    st.update(3, 4, 100)  # index 3 becomes 101
    assert st.query(0, 5) == 3


def test_max_explicit():
    st = LazySegmentTree([5, 3, 8, 1, 9], "max")
    assert st.query(0, 5) == 9
    st.update(0, 1, 100)  # index 0 becomes 105
    assert st.query(0, 5) == 105


def test_full_range_add():
    st = LazySegmentTree([1, 1, 1, 1], "sum")
    st.update(0, 4, 5)
    assert st.query(0, 4) == 24  # each becomes 6


def test_overlapping_updates():
    st = LazySegmentTree([0, 0, 0, 0, 0], "sum")
    st.update(0, 3, 1)
    st.update(2, 5, 2)
    # [1, 1, 3, 2, 2]
    assert st.query(0, 5) == 9
    assert st.query(2, 3) == 3


def test_empty_range_query():
    st = LazySegmentTree([1, 2, 3], "sum")
    assert st.query(2, 2) == 0
    assert LazySegmentTree([1, 2, 3], "min").query(1, 1) == float("inf")


def test_single_element():
    st = LazySegmentTree([7], "sum")
    st.update(0, 1, 3)
    assert st.query(0, 1) == 10


def test_len():
    assert len(LazySegmentTree([1, 2, 3, 4])) == 4


def test_bad_mode_raises():
    with pytest.raises(ValueError):
        LazySegmentTree([1], "product")


def test_out_of_bounds_raises():
    st = LazySegmentTree([1, 2, 3])
    with pytest.raises(IndexError):
        st.query(0, 5)
    with pytest.raises(IndexError):
        st.update(-1, 2, 1)


def test_non_power_of_two_length():
    # internal padding must not leak into queries over a 5-element (non-2^k) array
    st = LazySegmentTree([2, 4, 6, 8, 10], "sum")
    assert st.query(0, 5) == 30
    st.update(0, 5, 1)
    assert st.query(0, 5) == 35
    assert st.query(4, 5) == 11
