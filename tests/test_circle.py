"""Circle geometry: circumcircle, circle-line and circle-circle intersection."""

import math
import random

import pytest

from quantforge import (
    circle_from_3points,
    point_in_circle,
    circle_line_intersection,
    circle_circle_intersection,
)


def _on_circle(p, c, tol=1e-7):
    return abs(math.hypot(p[0] - c[0], p[1] - c[1]) - c[2]) < tol


def test_circumcircle():
    rng = random.Random(1)
    for _ in range(1000):
        pts = [(rng.uniform(-10, 10), rng.uniform(-10, 10)) for _ in range(3)]
        try:
            c = circle_from_3points(*pts)
        except ValueError:
            continue
        for p in pts:
            assert _on_circle(p, c)
    cx, cy, r = circle_from_3points((1, 0), (0, 1), (-1, 0))
    assert abs(cx) < 1e-9 and abs(cy) < 1e-9 and abs(r - 1) < 1e-9
    with pytest.raises(ValueError):
        circle_from_3points((0, 0), (1, 1), (2, 2))


def test_circle_line():
    c = (0, 0, 5)
    assert sorted(circle_line_intersection(c, (-10, 0), (10, 0))) == [(-5.0, 0.0), (5.0, 0.0)]
    assert circle_line_intersection(c, (-10, 5), (10, 5)) == [(0.0, 5.0)]      # tangent
    assert circle_line_intersection(c, (-10, 10), (10, 10)) == []             # miss

    rng = random.Random(2)
    for _ in range(2000):
        a = (rng.uniform(-10, 10), rng.uniform(-10, 10))
        b = (rng.uniform(-10, 10), rng.uniform(-10, 10))
        if a == b:
            continue
        for p in circle_line_intersection(c, a, b):
            assert _on_circle(p, c)


def test_circle_circle():
    pts = circle_circle_intersection((0, 0, 5), (8, 0, 5))
    assert len(pts) == 2 and all(_on_circle(p, (0, 0, 5)) and _on_circle(p, (8, 0, 5))
                                 for p in pts)
    assert circle_circle_intersection((0, 0, 1), (10, 0, 1)) == []       # separate
    assert circle_circle_intersection((0, 0, 5), (10, 0, 5)) == [(5.0, 0.0)]  # tangent
    with pytest.raises(ValueError):
        circle_circle_intersection((0, 0, 3), (0, 0, 3))                 # coincident


def test_point_in_circle():
    assert point_in_circle((1, 1), (0, 0, 2))
    assert not point_in_circle((3, 3), (0, 0, 2))
