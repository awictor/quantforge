"""Fenwick (binary indexed) tree and segment tree."""

import random

import pytest

from quantforge import FenwickTree, SegmentTree


def test_fenwick_vs_brute():
    rng = random.Random(1)
    for _ in range(200):
        n = rng.randint(1, 30)
        arr = [rng.uniform(-10, 10) for _ in range(n)]
        ft = FenwickTree(arr)
        for _ in range(20):
            i = rng.randrange(n)
            d = rng.uniform(-5, 5)
            arr[i] += d
            ft.update(i, d)
        for hi in range(n):
            assert abs(ft.prefix_sum(hi) - sum(arr[:hi + 1])) < 1e-6
        for _ in range(10):
            lo = rng.randrange(n)
            hi = rng.randint(lo, n - 1)
            assert abs(ft.range_sum(lo, hi) - sum(arr[lo:hi + 1])) < 1e-6


def test_fenwick_prefix_minus_one():
    assert FenwickTree([1, 2, 3]).prefix_sum(-1) == 0.0
    assert FenwickTree([1, 2, 3, 4, 5]).range_sum(1, 3) == 9.0


def test_segment_tree_sum():
    rng = random.Random(2)
    for _ in range(200):
        n = rng.randint(1, 30)
        arr = [rng.uniform(-10, 10) for _ in range(n)]
        st = SegmentTree(arr)
        for _ in range(15):
            i = rng.randrange(n)
            v = rng.uniform(-10, 10)
            arr[i] = v
            st.update(i, v)
        for _ in range(10):
            lo = rng.randrange(n)
            hi = rng.randint(lo, n - 1)
            assert abs(st.query(lo, hi) - sum(arr[lo:hi + 1])) < 1e-6


def test_segment_tree_min_max():
    rng = random.Random(3)
    for _ in range(200):
        n = rng.randint(1, 30)
        arr = [rng.randint(-100, 100) for _ in range(n)]
        smin = SegmentTree(arr, combine=min, identity=float("inf"))
        smax = SegmentTree(arr, combine=max, identity=float("-inf"))
        for _ in range(15):
            i = rng.randrange(n)
            v = rng.randint(-100, 100)
            arr[i] = v
            smin.update(i, v)
            smax.update(i, v)
        for _ in range(10):
            lo = rng.randrange(n)
            hi = rng.randint(lo, n - 1)
            assert smin.query(lo, hi) == min(arr[lo:hi + 1])
            assert smax.query(lo, hi) == max(arr[lo:hi + 1])


def test_validation():
    with pytest.raises(IndexError):
        FenwickTree(5).update(9, 1)
    with pytest.raises(ValueError):
        SegmentTree([])
    with pytest.raises(ValueError):
        FenwickTree([1, 2]).range_sum(1, 0)
