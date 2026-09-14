"""Polyline simplification (Ramer-Douglas-Peucker)."""

import random

import pytest

from quantforge import douglas_peucker
from quantforge.simplify import _perp_distance


def test_straight_line_collapses():
    line = [(i, 0) for i in range(10)]
    assert len(douglas_peucker(line, 0.01)) == 2


def test_invariant_and_subsequence():
    rng = random.Random(1)
    for _ in range(500):
        n = rng.randint(2, 40)
        pts = [(i, rng.uniform(-5, 5)) for i in range(n)]     # x == index
        eps = rng.uniform(0, 3)
        s = douglas_peucker(pts, eps)
        assert s[0] == pts[0] and s[-1] == pts[-1]
        for k in range(len(s) - 1):
            a, b = s[k], s[k + 1]
            for xi in range(int(a[0]) + 1, int(b[0])):
                assert _perp_distance(pts[xi], a, b) <= eps + 1e-9
        it = iter(pts)
        assert all(p in it for p in s)


def test_monotone_in_epsilon():
    zig = [(i, (-1) ** i) for i in range(20)]
    counts = [len(douglas_peucker(zig, e)) for e in (0, 0.5, 1.5, 3)]
    assert counts == sorted(counts, reverse=True)


def test_bump_keeps_peak():
    bump = [(0, 0), (1, 0), (2, 5), (3, 0), (4, 0)]
    assert (2, 5) in douglas_peucker(bump, 1)


def test_validation():
    with pytest.raises(ValueError):
        douglas_peucker([(0, 0), (1, 1)], -1)
