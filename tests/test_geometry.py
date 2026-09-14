"""Planar computational geometry."""

import math
import random

import pytest

from quantforge import (
    convex_hull,
    polygon_area,
    polygon_centroid,
    point_in_polygon,
    closest_pair,
)


def test_convex_hull_of_square_with_interior():
    pts = [(0, 0), (1, 0), (1, 1), (0, 1), (0.5, 0.5), (0.3, 0.7), (0.9, 0.1)]
    h = convex_hull(pts)
    assert set(h) == {(0, 0), (1, 0), (1, 1), (0, 1)}


def test_all_points_inside_their_hull():
    rng = random.Random(1)
    cloud = [(rng.uniform(0, 10), rng.uniform(0, 10)) for _ in range(200)]
    h = convex_hull(cloud)
    assert all(point_in_polygon(p, h) for p in cloud)


def test_polygon_area():
    assert abs(polygon_area([(0, 0), (1, 0), (1, 1), (0, 1)]) - 1.0) < 1e-12
    assert abs(polygon_area([(0, 0), (4, 0), (0, 3)]) - 6.0) < 1e-12
    # orientation independent
    assert abs(polygon_area([(0, 0), (0, 1), (1, 1), (1, 0)]) - 1.0) < 1e-12


def test_polygon_centroid():
    cx, cy = polygon_centroid([(0, 0), (2, 0), (2, 2), (0, 2)])
    assert abs(cx - 1.0) < 1e-12 and abs(cy - 1.0) < 1e-12
    cx, cy = polygon_centroid([(0, 0), (6, 0), (0, 3)])
    assert abs(cx - 2.0) < 1e-12 and abs(cy - 1.0) < 1e-12


def test_point_in_polygon():
    sq = [(0, 0), (4, 0), (4, 4), (0, 4)]
    assert point_in_polygon((2, 2), sq)
    assert not point_in_polygon((5, 5), sq)
    assert point_in_polygon((0, 2), sq)          # on edge
    # concave L-shape
    L = [(0, 0), (4, 0), (4, 2), (2, 2), (2, 4), (0, 4)]
    assert point_in_polygon((1, 3), L)
    assert not point_in_polygon((3, 3), L)       # in the notch


def test_closest_pair_matches_brute_force():
    rng = random.Random(2)

    def brute(pts):
        best = None
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                d = math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1])
                if best is None or d < best:
                    best = d
        return best

    for _ in range(50):
        cp = [(rng.uniform(0, 100), rng.uniform(0, 100)) for _ in range(60)]
        _, _, d = closest_pair(cp)
        assert abs(d - brute(cp)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        polygon_area([(0, 0), (1, 1)])           # fewer than 3 vertices
    with pytest.raises(ValueError):
        point_in_polygon((0, 0), [(0, 0), (1, 1)])
    with pytest.raises(ValueError):
        closest_pair([(0, 0)])                   # fewer than 2 points
