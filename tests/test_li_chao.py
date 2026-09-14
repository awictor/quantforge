"""Tests for the Li Chao tree, cross-checked against a brute min/max over all lines."""

import random

import pytest

from quantforge.li_chao import LiChaoTree


def test_fuzz_min_and_max_envelopes():
    rng = random.Random(321)
    for maximize in (False, True):
        pick = max if maximize else min
        for _ in range(3000):
            xs = sorted(set(rng.randint(-50, 50) for _ in range(rng.randint(1, 20))))
            t = LiChaoTree(xs, maximize=maximize)
            lines = []
            for _ in range(rng.randint(1, 15)):
                m = rng.randint(-10, 10)
                b = rng.randint(-50, 50)
                t.add_line(m, b)
                lines.append((m, b))
                for x in xs:
                    assert t.query(x) == pick(m2 * x + b2 for m2, b2 in lines)


def test_min_envelope_explicit():
    t = LiChaoTree([0, 1, 2, 3, 4])
    t.add_line(0, 5)     # y = 5
    t.add_line(2, -1)    # y = 2x - 1
    t.add_line(-2, 8)    # y = -2x + 8
    assert t.query(0) == -1  # min(5, -1, 8)
    assert t.query(4) == 0   # min(5, 7, 0)


def test_max_envelope_explicit():
    t = LiChaoTree([0, 1, 2, 3, 4], maximize=True)
    t.add_line(0, 5)
    t.add_line(2, -1)
    t.add_line(-2, 8)
    assert t.query(0) == 8   # max(5, -1, 8)
    assert t.query(4) == 7   # max(5, 7, 0)


def test_single_line():
    t = LiChaoTree([1, 2, 3])
    t.add_line(3, 1)  # y = 3x + 1
    assert t.query(2) == 7


def test_single_point():
    t = LiChaoTree([5])
    t.add_line(1, 1)
    t.add_line(2, -3)
    assert t.query(5) == min(6, 7)


def test_parallel_lines():
    t = LiChaoTree([0, 1, 2])
    t.add_line(1, 10)
    t.add_line(1, 3)  # same slope, lower intercept dominates
    assert t.query(0) == 3
    assert t.query(2) == 5


def test_insertion_order_independent():
    xs = [0, 1, 2, 3, 4, 5]
    lines = [(0, 5), (2, -1), (-2, 8), (1, 0), (-1, 6)]
    a = LiChaoTree(xs)
    for m, b in lines:
        a.add_line(m, b)
    b_tree = LiChaoTree(xs)
    for m, b in reversed(lines):
        b_tree.add_line(m, b)
    for x in xs:
        assert a.query(x) == b_tree.query(x)


def test_dp_optimization_pattern():
    # dp[i] = min_j (m_j * x_i + b_j): the classic use case
    xs = [1, 2, 3, 4, 5]
    t = LiChaoTree(xs)
    t.add_line(3, 0)
    t.add_line(1, 4)
    t.add_line(-1, 10)
    assert t.query(5) == min(15, 9, 5)


def test_empty_construction_raises():
    with pytest.raises(ValueError):
        LiChaoTree([])


def test_query_unknown_x_raises():
    t = LiChaoTree([1, 2, 3])
    t.add_line(1, 0)
    with pytest.raises(ValueError):
        t.query(9)
