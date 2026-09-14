"""Planar geometry III: diameter, bounding box, minimum enclosing circle."""

import math
import random

import pytest

from quantforge import bounding_box, polygon_diameter, min_enclosing_circle


def test_bounding_box():
    assert bounding_box([(1, 2), (3, -1), (0, 5), (4, 4)]) == (0, -1, 4, 5)
    assert bounding_box([(2, 2)]) == (2, 2, 2, 2)


def test_diameter_matches_brute_force():
    rng = random.Random(1)

    def brute(pts):
        best = 0.0
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                best = max(best, math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1]))
        return best

    for _ in range(50):
        pts = [(rng.uniform(0, 100), rng.uniform(0, 100)) for _ in range(40)]
        _, _, d = polygon_diameter(pts)
        assert abs(d - brute(pts)) < 1e-9


def test_mec_contains_all_and_minimal():
    rng = random.Random(2)

    def contains(pts, c):
        return all(math.hypot(p[0] - c[0], p[1] - c[1]) <= c[2] + 1e-7 for p in pts)

    for _ in range(50):
        pts = [(rng.uniform(0, 100), rng.uniform(0, 100)) for _ in range(30)]
        c = min_enclosing_circle(pts)
        assert contains(pts, c)
        assert not contains(pts, (c[0], c[1], c[2] * 0.99))    # can't shrink


def test_mec_known_cases():
    c = min_enclosing_circle([(0, 0), (2, 0)])
    assert abs(c[0] - 1) < 1e-9 and abs(c[1]) < 1e-9 and abs(c[2] - 1) < 1e-9
    c = min_enclosing_circle([(0, 0), (2, 0), (2, 2), (0, 2)])
    assert abs(c[0] - 1) < 1e-9 and abs(c[1] - 1) < 1e-9 and abs(c[2] - math.sqrt(2)) < 1e-9
    assert min_enclosing_circle([(5, 5)]) == (5, 5, 0.0)


def test_validation():
    with pytest.raises(ValueError):
        bounding_box([])
    with pytest.raises(ValueError):
        polygon_diameter([(0, 0)])
    with pytest.raises(ValueError):
        min_enclosing_circle([])
