"""Tests for SparseTable O(1) idempotent range queries, cross-checked against brute force."""

import math
import random

import pytest

from quantforge.sparse_table import (
    SparseTable,
    range_min_query,
    range_max_query,
    range_gcd_query,
)


def test_fuzz_min_max_vs_brute():
    rng = random.Random(41)
    for _ in range(3000):
        n = rng.randint(1, 60)
        a = [rng.randint(-100, 100) for _ in range(n)]
        st_min = range_min_query(a)
        st_max = range_max_query(a)
        for _ in range(5):
            lo = rng.randint(0, n - 1)
            hi = rng.randint(lo, n - 1)
            assert st_min.query(lo, hi) == min(a[lo:hi + 1])
            assert st_max.query(lo, hi) == max(a[lo:hi + 1])


def test_fuzz_gcd_vs_brute():
    rng = random.Random(42)
    for _ in range(2000):
        n = rng.randint(1, 50)
        a = [rng.randint(1, 60) for _ in range(n)]
        st = range_gcd_query(a)
        for _ in range(5):
            lo = rng.randint(0, n - 1)
            hi = rng.randint(lo, n - 1)
            g = 0
            for x in a[lo:hi + 1]:
                g = math.gcd(g, x)
            assert st.query(lo, hi) == g


def test_single_element():
    st = range_min_query([5])
    assert st.query(0, 0) == 5


def test_full_and_sub_ranges():
    st = range_max_query([3, 1, 4, 1, 5, 9, 2, 6])
    assert st.query(0, 7) == 9
    assert st.query(1, 3) == 4
    assert st.query(6, 7) == 6


def test_min_known_values():
    st = range_min_query([3, 1, 4, 1, 5, 9, 2, 6])
    assert st.query(0, 2) == 1
    assert st.query(4, 7) == 2
    assert st.query(2, 2) == 4


def test_custom_idempotent_combine_bitwise_or():
    st = SparseTable([1, 2, 4, 8], combine=lambda a, b: a | b)
    assert st.query(0, 3) == 15
    assert st.query(1, 2) == 6


def test_len():
    assert len(range_min_query([1, 2, 3])) == 3
    assert len(SparseTable([])) == 0


def test_empty_query_raises():
    st = range_min_query([1, 2])
    with pytest.raises(IndexError):
        st.query(0, 5)
    with pytest.raises(IndexError):
        st.query(-1, 1)
    with pytest.raises(IndexError):
        st.query(1, 0)


def test_power_of_two_boundaries():
    # exercise block lengths that are and are not powers of two
    a = list(range(16, 0, -1))  # 16, 15, ..., 1
    st = range_min_query(a)
    for lo in range(len(a)):
        for hi in range(lo, len(a)):
            assert st.query(lo, hi) == min(a[lo:hi + 1])


def test_query_o1_after_build_matches_recompute():
    a = [5, 3, 8, 8, 3, 1, 7, 2, 9, 4]
    st_min = range_min_query(a)
    st_max = range_max_query(a)
    for lo in range(len(a)):
        for hi in range(lo, len(a)):
            assert st_min.query(lo, hi) == min(a[lo:hi + 1])
            assert st_max.query(lo, hi) == max(a[lo:hi + 1])
