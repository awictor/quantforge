"""Tests for the minimum-area enclosing rectangle (rotating calipers)."""

import math
import random

import pytest

from quantforge.min_area_rect import min_area_rectangle


def _close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _brute_min_area(points):
    best = float("inf")
    for k in range(360):
        th = k * math.pi / 360
        c, s = math.cos(th), math.sin(th)
        us = [p[0] * c + p[1] * s for p in points]
        vs = [-p[0] * s + p[1] * c for p in points]
        best = min(best, (max(us) - min(us)) * (max(vs) - min(vs)))
    return best


def _point_in_rect(p, corners, tol=1e-6):
    n = len(corners)
    sign = None
    for i in range(n):
        ax, ay = corners[i]
        bx, by = corners[(i + 1) % n]
        cr = (bx - ax) * (p[1] - ay) - (by - ay) * (p[0] - ax)
        if abs(cr) < tol:
            continue
        s = cr > 0
        if sign is None:
            sign = s
        elif s != sign:
            return False
    return True


def test_fuzz_optimal_and_contains():
    rng = random.Random(491)
    for _ in range(2000):
        n = rng.randint(3, 20)
        pts = [(rng.uniform(-10, 10), rng.uniform(-10, 10)) for _ in range(n)]
        r = min_area_rectangle(pts)
        assert r["area"] <= _brute_min_area(pts) + 1e-6
        for p in pts:
            assert _point_in_rect(p, r["corners"], 1e-6)
        assert r["width"] >= r["height"] - 1e-9


def test_axis_aligned_square():
    r = min_area_rectangle([(0, 0), (2, 0), (2, 2), (0, 2)])
    assert _close(r["area"], 4.0)
    assert _close(r["width"], 2.0)
    assert _close(r["height"], 2.0)


def test_rotated_diamond():
    r = min_area_rectangle([(1, 0), (2, 1), (1, 2), (0, 1)])
    assert _close(r["area"], 2.0)


def test_rectangle():
    r = min_area_rectangle([(0, 0), (4, 0), (4, 1), (0, 1)])
    assert _close(r["area"], 4.0)
    assert _close(r["width"], 4.0)
    assert _close(r["height"], 1.0)


def test_interior_points_ignored():
    # interior points do not change the bounding rectangle
    r = min_area_rectangle([(0, 0), (4, 0), (4, 3), (0, 3), (1, 1), (2, 2), (3, 1)])
    assert _close(r["area"], 12.0)


def test_collinear_zero_area():
    r = min_area_rectangle([(0, 0), (1, 1), (2, 2)])
    assert _close(r["area"], 0.0)


def test_single_point():
    r = min_area_rectangle([(3, 3)])
    assert r["area"] == 0.0
    assert r["corners"] == [(3, 3)] * 4


def test_two_points():
    r = min_area_rectangle([(0, 0), (3, 4)])
    assert r["area"] == 0.0
    assert _close(r["width"], 5.0)


def test_corners_form_rectangle():
    r = min_area_rectangle([(0, 0), (5, 0), (5, 2), (0, 2), (2, 1)])
    c = r["corners"]
    # opposite sides equal length
    d01 = math.hypot(c[1][0] - c[0][0], c[1][1] - c[0][1])
    d23 = math.hypot(c[3][0] - c[2][0], c[3][1] - c[2][1])
    assert _close(d01, d23)


def test_empty_raises():
    with pytest.raises(ValueError):
        min_area_rectangle([])


def test_rotated_rectangle_recovered():
    # a 6x2 rectangle rotated 30 deg: min area must still be 12
    th = math.pi / 6
    c, s = math.cos(th), math.sin(th)
    base = [(0, 0), (6, 0), (6, 2), (0, 2)]
    rot = [(x * c - y * s, x * s + y * c) for x, y in base]
    r = min_area_rectangle(rot)
    assert _close(r["area"], 12.0, 1e-6)
