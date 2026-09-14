"""Tests for prefix-sum structures, cross-checked against brute-force sums."""

import random

import pytest

from quantforge.prefix_sum import PrefixSum1D, PrefixSum2D, DifferenceArray


def test_fuzz_1d_range_sum():
    rng = random.Random(171)
    for _ in range(3000):
        n = rng.randint(0, 30)
        a = [rng.randint(-20, 20) for _ in range(n)]
        ps = PrefixSum1D(a)
        for _ in range(5):
            lo = rng.randint(0, n)
            hi = rng.randint(lo, n)
            assert ps.range_sum(lo, hi) == sum(a[lo:hi])
        assert ps.total() == sum(a)


def test_fuzz_2d_rectangle_sum():
    rng = random.Random(172)
    for _ in range(2000):
        R = rng.randint(0, 8)
        C = rng.randint(0, 8) if R else 0
        grid = [[rng.randint(-10, 10) for _ in range(C)] for _ in range(R)]
        ps = PrefixSum2D(grid)
        for _ in range(5):
            r0 = rng.randint(0, R)
            r1 = rng.randint(r0, R)
            c0 = rng.randint(0, C)
            c1 = rng.randint(c0, C)
            brute = sum(grid[i][j] for i in range(r0, r1) for j in range(c0, c1))
            assert ps.range_sum(r0, c0, r1, c1) == brute


def test_fuzz_difference_array():
    rng = random.Random(173)
    for _ in range(3000):
        n = rng.randint(0, 25)
        base = [rng.randint(-5, 5) for _ in range(n)]
        da = DifferenceArray(base)
        ref = base[:]
        for _ in range(rng.randint(0, 10)):
            lo = rng.randint(0, n)
            hi = rng.randint(lo, n)
            d = rng.randint(-7, 7)
            da.add(lo, hi, d)
            for i in range(lo, hi):
                ref[i] += d
        assert da.result() == ref


def test_1d_explicit():
    ps = PrefixSum1D([1, 2, 3, 4, 5])
    assert ps.range_sum(1, 4) == 9  # 2 + 3 + 4
    assert ps.total() == 15
    assert ps.prefix(2) == 3
    assert ps.range_sum(0, 0) == 0
    assert len(ps) == 5


def test_2d_explicit():
    grid = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    ps = PrefixSum2D(grid)
    assert ps.range_sum(0, 0, 3, 3) == 45
    assert ps.range_sum(1, 1, 3, 3) == 28  # 5+6+8+9
    assert ps.range_sum(0, 0, 1, 1) == 1
    assert ps.shape == (3, 3)
    assert ps.total() == 45


def test_difference_array_explicit():
    da = DifferenceArray(5)
    da.add(1, 4, 10)
    da.add(0, 2, 1)
    assert da.result() == [1, 11, 10, 10, 0]


def test_difference_array_from_int_length():
    da = DifferenceArray(3)
    assert da.result() == [0, 0, 0]
    da.add(0, 3, 5)
    assert da.result() == [5, 5, 5]


def test_empty_structures():
    assert PrefixSum1D([]).total() == 0
    assert PrefixSum2D([]).total() == 0
    assert DifferenceArray(0).result() == []
    assert DifferenceArray([]).result() == []


def test_1d_out_of_bounds_raises():
    ps = PrefixSum1D([1, 2, 3])
    with pytest.raises(IndexError):
        ps.range_sum(0, 5)
    with pytest.raises(IndexError):
        ps.range_sum(2, 1)
    with pytest.raises(IndexError):
        ps.prefix(9)


def test_2d_ragged_raises():
    with pytest.raises(ValueError):
        PrefixSum2D([[1, 2], [3]])


def test_2d_out_of_bounds_raises():
    ps = PrefixSum2D([[1, 2], [3, 4]])
    with pytest.raises(IndexError):
        ps.range_sum(0, 0, 3, 3)


def test_difference_array_out_of_bounds_raises():
    da = DifferenceArray(3)
    with pytest.raises(IndexError):
        da.add(0, 5, 1)
    with pytest.raises(IndexError):
        da.add(2, 1, 1)


def test_difference_array_negative_length_raises():
    with pytest.raises(ValueError):
        DifferenceArray(-1)


def test_add_returns_self():
    da = DifferenceArray(5)
    assert da.add(0, 2, 1) is da


def test_single_row_2d():
    ps = PrefixSum2D([[1, 2, 3, 4]])
    assert ps.range_sum(0, 1, 1, 3) == 5  # cols 1..2 of the single row
