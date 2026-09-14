"""Tests for the Voronoi diagram (Delaunay dual), cross-checked against brute nearest-site."""

import math
import random

import pytest

from quantforge.voronoi import voronoi_vertices, delaunay_neighbors, nearest_site
from quantforge.delaunay import delaunay_triangulation


def _close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _random_pts(rng, n):
    return list(dict.fromkeys(
        (round(rng.uniform(0, 10), 3), round(rng.uniform(0, 10), 3)) for _ in range(n)
    ))


def test_fuzz_vertices_equidistant_and_empty():
    rng = random.Random(511)
    for _ in range(1500):
        pts = _random_pts(rng, rng.randint(3, 15))
        if len(pts) < 3:
            continue
        try:
            tris = delaunay_triangulation(pts)
            verts = voronoi_vertices(pts)
        except ValueError:
            continue
        for (i, j, k), v in zip(tris, verts):
            d = [math.hypot(v[0] - pts[m][0], v[1] - pts[m][1]) for m in (i, j, k)]
            assert _close(d[0], d[1], 1e-5) and _close(d[1], d[2], 1e-5)
            r = d[0]
            for m, p in enumerate(pts):
                if m in (i, j, k):
                    continue
                assert math.hypot(v[0] - p[0], v[1] - p[1]) >= r - 1e-5


def test_fuzz_adjacency_symmetric():
    rng = random.Random(513)
    for _ in range(1000):
        pts = _random_pts(rng, rng.randint(3, 15))
        if len(pts) < 3:
            continue
        try:
            adj = delaunay_neighbors(pts)
        except ValueError:
            continue
        for i in adj:
            for j in adj[i]:
                assert i in adj[j]


def test_fuzz_nearest_site_vs_brute():
    rng = random.Random(512)
    for _ in range(3000):
        pts = [(rng.uniform(0, 10), rng.uniform(0, 10)) for _ in range(rng.randint(1, 12))]
        q = (rng.uniform(0, 10), rng.uniform(0, 10))
        ns = nearest_site(pts, q)
        brute = min(range(len(pts)), key=lambda i: (pts[i][0] - q[0]) ** 2 + (pts[i][1] - q[1]) ** 2)
        d_ns = (pts[ns][0] - q[0]) ** 2 + (pts[ns][1] - q[1]) ** 2
        d_br = (pts[brute][0] - q[0]) ** 2 + (pts[brute][1] - q[1]) ** 2
        assert d_ns <= d_br + 1e-12


def test_square_voronoi_center():
    pts = [(0, 0), (2, 0), (2, 2), (0, 2)]
    verts = voronoi_vertices(pts)
    for v in verts:
        assert _close(v[0], 1.0)
        assert _close(v[1], 1.0)


def test_square_adjacency():
    pts = [(0, 0), (2, 0), (2, 2), (0, 2)]
    adj = delaunay_neighbors(pts)
    # every corner is adjacent to the two edge-sharing corners
    assert 1 in adj[0] and 3 in adj[0]
    assert 0 in adj[1] and 2 in adj[1]


def test_nearest_site_exact_hit():
    pts = [(0, 0), (5, 5), (10, 0)]
    assert nearest_site(pts, (5, 5)) == 1
    assert nearest_site(pts, (0.1, 0.1)) == 0


def test_single_site():
    assert nearest_site([(3, 3)], (100, 100)) == 0


def test_triangle_one_vertex():
    verts = voronoi_vertices([(0, 0), (4, 0), (2, 3)])
    assert len(verts) == 1
    # equidistant from all three
    v = verts[0]
    d = [math.hypot(v[0] - p[0], v[1] - p[1]) for p in [(0, 0), (4, 0), (2, 3)]]
    assert _close(d[0], d[1]) and _close(d[1], d[2])


def test_empty_nearest_raises():
    with pytest.raises(ValueError):
        nearest_site([], (0, 0))


def test_too_few_points_raises():
    with pytest.raises(ValueError):
        voronoi_vertices([(0, 0), (1, 1)])
