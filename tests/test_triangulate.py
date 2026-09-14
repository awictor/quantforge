"""Polygon triangulation by ear clipping, orientation, convexity."""

import math
import random

import pytest

from quantforge import (
    ear_clipping_triangulate,
    signed_area,
    is_clockwise,
    is_convex_polygon,
    polygon_area,
)


def _tri_area(t):
    (ax, ay), (bx, by), (cx, cy) = t
    return abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay)) / 2


def test_square():
    tris = ear_clipping_triangulate([(0, 0), (4, 0), (4, 4), (0, 4)])
    assert len(tris) == 2
    assert abs(sum(_tri_area(t) for t in tris) - 16) < 1e-9


def test_concave_L_shape():
    L = [(0, 0), (4, 0), (4, 2), (2, 2), (2, 4), (0, 4)]
    tris = ear_clipping_triangulate(L)
    assert len(tris) == len(L) - 2
    assert abs(sum(_tri_area(t) for t in tris) - polygon_area(L)) < 1e-9


def test_random_convex_polygons():
    for seed in range(300):
        rng = random.Random(seed * 7 + 1)
        n = rng.randint(3, 12)
        angles = sorted(rng.uniform(0, 2 * math.pi) for _ in range(n))
        poly = [(2 * math.cos(a), 2 * math.sin(a)) for a in angles]
        if len(set(poly)) < 3:
            continue
        tris = ear_clipping_triangulate(poly)
        assert len(tris) == len(poly) - 2
        assert abs(sum(_tri_area(t) for t in tris) - polygon_area(poly)) < 1e-6


def test_orientation_and_convexity():
    sq = [(0, 0), (4, 0), (4, 4), (0, 4)]
    assert signed_area(sq) == 16.0
    assert not is_clockwise(sq)
    assert is_convex_polygon(sq)
    assert is_clockwise([(0, 0), (0, 4), (4, 4), (4, 0)])
    assert not is_convex_polygon([(0, 0), (4, 0), (4, 2), (2, 2), (2, 4), (0, 4)])


def test_validation():
    with pytest.raises(ValueError):
        ear_clipping_triangulate([(0, 0), (1, 1)])
