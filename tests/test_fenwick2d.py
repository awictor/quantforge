"""Tests for the 2-D Fenwick tree, cross-checked against a brute grid."""

import random

import pytest

from quantforge.fenwick2d import FenwickTree2D


def _brute(g, r0, c0, r1, c1):
    return sum(g[i][j] for i in range(r0, r1) for j in range(c0, c1))


def test_fuzz_interleaved_updates():
    rng = random.Random(331)
    for _ in range(3000):
        R = rng.randint(1, 8)
        C = rng.randint(1, 8)
        g = [[0] * C for _ in range(R)]
        ft = FenwickTree2D(R, C)
        for _ in range(15):
            if rng.random() < 0.5:
                r = rng.randint(0, R - 1)
                c = rng.randint(0, C - 1)
                d = rng.randint(-10, 10)
                ft.add(r, c, d)
                g[r][c] += d
            else:
                r0 = rng.randint(0, R)
                r1 = rng.randint(r0, R)
                c0 = rng.randint(0, C)
                c1 = rng.randint(c0, C)
                assert ft.range_sum(r0, c0, r1, c1) == _brute(g, r0, c0, r1, c1)


def test_fuzz_init_from_grid():
    rng = random.Random(332)
    for _ in range(1000):
        R = rng.randint(1, 6)
        C = rng.randint(1, 6)
        g = [[rng.randint(-5, 5) for _ in range(C)] for _ in range(R)]
        ft = FenwickTree2D(grid=g)
        assert ft.shape == (R, C)
        for _ in range(5):
            r0 = rng.randint(0, R)
            r1 = rng.randint(r0, R)
            c0 = rng.randint(0, C)
            c1 = rng.randint(c0, C)
            assert ft.range_sum(r0, c0, r1, c1) == _brute(g, r0, c0, r1, c1)


def test_explicit():
    ft = FenwickTree2D(3, 3)
    for i in range(3):
        for j in range(3):
            ft.add(i, j, i * 3 + j + 1)  # values 1..9
    assert ft.range_sum(0, 0, 3, 3) == 45
    assert ft.range_sum(1, 1, 3, 3) == 28  # 5+6+8+9
    assert ft.prefix_sum(2, 2) == 12       # 1+2+4+5


def test_update_then_query():
    ft = FenwickTree2D(3, 3)
    ft.add(1, 1, 100)
    assert ft.range_sum(1, 1, 2, 2) == 100
    assert ft.range_sum(0, 0, 3, 3) == 100


def test_empty_range():
    ft = FenwickTree2D(3, 3)
    ft.add(0, 0, 5)
    assert ft.range_sum(1, 1, 1, 1) == 0
    assert ft.prefix_sum(0, 0) == 0


def test_single_cell_grid():
    ft = FenwickTree2D(1, 1)
    ft.add(0, 0, 42)
    assert ft.range_sum(0, 0, 1, 1) == 42


def test_negative_deltas():
    ft = FenwickTree2D(2, 2)
    ft.add(0, 0, 10)
    ft.add(0, 0, -3)
    assert ft.range_sum(0, 0, 1, 1) == 7


def test_shape():
    assert FenwickTree2D(4, 7).shape == (4, 7)
    assert FenwickTree2D(grid=[[1, 2, 3]]).shape == (1, 3)


def test_out_of_bounds_add_raises():
    ft = FenwickTree2D(2, 2)
    with pytest.raises(IndexError):
        ft.add(5, 0, 1)
    with pytest.raises(IndexError):
        ft.add(0, -1, 1)


def test_out_of_bounds_query_raises():
    ft = FenwickTree2D(2, 2)
    with pytest.raises(IndexError):
        ft.range_sum(0, 0, 5, 5)


def test_ragged_grid_raises():
    with pytest.raises(ValueError):
        FenwickTree2D(grid=[[1, 2], [3]])


def test_missing_dimensions_raises():
    with pytest.raises(ValueError):
        FenwickTree2D(3)  # no cols, no grid
