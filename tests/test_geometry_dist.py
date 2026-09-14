"""Point-to-line and point-to-segment distances."""

import math
import random

import pytest

from quantforge import (
    point_to_line_distance,
    closest_point_on_segment,
    point_segment_distance,
    point_polyline_distance,
)


def test_perpendicular_distance():
    assert point_to_line_distance((0, 3), (0, 0), (1, 0)) == 3
    assert abs(point_to_line_distance((1, 1), (0, 0), (2, 2))) < 1e-12
    assert abs(point_to_line_distance((0, 0), (2, 0), (0, 2)) - math.sqrt(2)) < 1e-9


def test_segment_clamping():
    assert point_segment_distance((5, 0), (0, 0), (3, 0)) == 2      # beyond endpoint
    assert point_segment_distance((1, 4), (0, 0), (3, 0)) == 4      # perpendicular
    assert point_segment_distance((1.5, 0), (0, 0), (3, 0)) == 0    # on segment


def test_closest_point_vs_brute():
    rng = random.Random(1)
    for _ in range(2000):
        a = (rng.uniform(-10, 10), rng.uniform(-10, 10))
        b = (rng.uniform(-10, 10), rng.uniform(-10, 10))
        if a == b:
            continue
        p = (rng.uniform(-10, 10), rng.uniform(-10, 10))
        cp = closest_point_on_segment(p, a, b)
        d = math.hypot(p[0] - cp[0], p[1] - cp[1])
        brute = min(math.hypot(p[0] - (a[0] + t / 1000 * (b[0] - a[0])),
                               p[1] - (a[1] + t / 1000 * (b[1] - a[1])))
                    for t in range(1001))
        assert d <= brute + 1e-9
        assert abs(point_segment_distance(p, a, b) - d) < 1e-12


def test_polyline():
    poly = [(0, 0), (2, 0), (2, 2), (4, 2)]
    assert abs(point_polyline_distance((3, 1), poly) - 1) < 1e-9
    assert abs(point_polyline_distance((1, -1), poly) - 1) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        point_to_line_distance((0, 0), (1, 1), (1, 1))
    with pytest.raises(ValueError):
        point_polyline_distance((0, 0), [(1, 1)])
