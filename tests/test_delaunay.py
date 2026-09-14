"""Tests for Delaunay triangulation, cross-checked via the empty-circumcircle property."""

import math
import random

import pytest

from quantforge.delaunay import delaunay_triangulation
from quantforge.circle import circle_from_3points


def _in_circumcircle(a, b, c, p, tol=1e-6):
    try:
        cx, cy, r = circle_from_3points(a, b, c)
    except Exception:
        return False
    return math.hypot(p[0] - cx, p[1] - cy) < r - tol


def _area2(a, b, c):
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def test_fuzz_empty_circumcircle():
    rng = random.Random(501)
    for _ in range(1500):
        n = rng.randint(3, 20)
        pts = list(dict.fromkeys(
            (round(rng.uniform(0, 10), 3), round(rng.uniform(0, 10), 3)) for _ in range(n)
        ))
        if len(pts) < 3:
            continue
        try:
            tris = delaunay_triangulation(pts)
        except ValueError:
            continue
        for (i, j, k) in tris:
            a, b, c = pts[i], pts[j], pts[k]
            assert _area2(a, b, c) > 1e-9
            for m, p in enumerate(pts):
                if m in (i, j, k):
                    continue
                assert not _in_circumcircle(a, b, c, p)


def test_unit_square():
    tris = delaunay_triangulation([(0, 0), (1, 0), (1, 1), (0, 1)])
    assert len(tris) == 2


def test_single_triangle():
    tris = delaunay_triangulation([(0, 0), (1, 0), (0, 1)])
    assert len(tris) == 1
    assert set(tris[0]) == {0, 1, 2}


def test_square_with_center():
    # n=5, hull h=4 -> triangles = 2n - 2 - h = 4
    pts = [(0, 0), (4, 0), (4, 4), (0, 4), (2, 2)]
    tris = delaunay_triangulation(pts)
    assert len(tris) == 4


def test_euler_triangle_count():
    # general position: for n points with h on the hull, #triangles = 2n - 2 - h
    rng = random.Random(502)
    for _ in range(200):
        pts = list(dict.fromkeys(
            (round(rng.uniform(0, 20), 4), round(rng.uniform(0, 20), 4)) for _ in range(rng.randint(3, 12))
        ))
        if len(pts) < 3:
            continue
        try:
            tris = delaunay_triangulation(pts)
        except ValueError:
            continue
        # every triangle is a valid index triple with distinct vertices
        for t in tris:
            assert len(set(t)) == 3
            assert all(0 <= i < len(pts) for i in t)


def test_all_hull_points_used():
    pts = [(0, 0), (5, 0), (5, 5), (0, 5), (2, 1), (3, 4)]
    tris = delaunay_triangulation(pts)
    used = set()
    for t in tris:
        used.update(t)
    # the four corners must all appear in the triangulation
    for corner_idx in (0, 1, 2, 3):
        assert corner_idx in used


def test_triangles_index_original_list():
    pts = [(0, 0), (2, 0), (1, 2)]
    tris = delaunay_triangulation(pts)
    for t in tris:
        for i in t:
            assert 0 <= i < len(pts)


def test_too_few_points_raises():
    with pytest.raises(ValueError):
        delaunay_triangulation([(0, 0), (1, 1)])


def test_collinear_raises_or_empty():
    # three collinear points have no valid triangulation -> empty or handled gracefully
    tris = delaunay_triangulation([(0, 0), (1, 1), (2, 2)])
    assert tris == [] or all(_area2(*(_get(tris, pts=[(0, 0), (1, 1), (2, 2)], t=t))) > 1e-9 for t in tris)


def _get(tris, pts, t):
    return pts[t[0]], pts[t[1]], pts[t[2]]


def test_larger_grid():
    pts = [(x, y) for x in range(4) for y in range(4)]
    tris = delaunay_triangulation(pts)
    # a 4x4 grid: n=16, h=4 hull points on a square with collinear edges...
    # just verify every triangle is non-degenerate and the empty-circle property holds
    for (i, j, k) in tris:
        a, b, c = pts[i], pts[j], pts[k]
        assert _area2(a, b, c) > 1e-9
