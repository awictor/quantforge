"""Planar geometry II: segment intersection and convex polygon clipping."""

import pytest

from quantforge import (
    segments_intersect,
    segment_intersection,
    polygon_perimeter,
    clip_polygon,
    polygon_area,
)


def test_crossing_segments():
    assert segments_intersect((0, 0), (1, 1), (0, 1), (1, 0))
    pt = segment_intersection((0, 0), (1, 1), (0, 1), (1, 0))
    assert abs(pt[0] - 0.5) < 1e-12 and abs(pt[1] - 0.5) < 1e-12


def test_parallel_and_disjoint():
    assert segment_intersection((0, 0), (1, 0), (0, 1), (1, 1)) is None
    assert not segments_intersect((0, 0), (1, 0), (2, 2), (3, 3))


def test_endpoint_and_t_junction():
    assert segments_intersect((0, 0), (1, 1), (1, 1), (2, 0))       # shared endpoint
    assert segments_intersect((0, 0), (2, 0), (1, 0), (1, 5))       # T-junction


def test_perimeter():
    assert abs(polygon_perimeter([(0, 0), (1, 0), (1, 1), (0, 1)]) - 4.0) < 1e-12
    assert abs(polygon_perimeter([(0, 0), (4, 0), (0, 3)]) - 12.0) < 1e-12   # 3-4-5


def test_clip_square_by_window():
    sq = [(0, 0), (4, 0), (4, 4), (0, 4)]
    win = [(1, 1), (3, 1), (3, 3), (1, 3)]
    clipped = clip_polygon(sq, win)
    assert abs(polygon_area(clipped) - 4.0) < 1e-9           # the 2x2 window


def test_clip_contained_and_outside():
    tri = [(1, 1), (3, 1), (2, 2.5)]
    inside = clip_polygon(tri, [(0, 0), (5, 0), (5, 5), (0, 5)])
    assert abs(polygon_area(inside) - polygon_area(tri)) < 1e-9
    assert clip_polygon([(10, 10), (11, 10), (11, 11)],
                        [(1, 1), (3, 1), (3, 3), (1, 3)]) == []


def test_clip_partial_overlap():
    sq = [(0, 0), (4, 0), (4, 4), (0, 4)]
    win = [(2, 2), (6, 2), (6, 6), (2, 6)]
    assert abs(polygon_area(clip_polygon(sq, win)) - 4.0) < 1e-9   # [2,4]^2


def test_validation():
    with pytest.raises(ValueError):
        polygon_perimeter([(0, 0)])
    with pytest.raises(ValueError):
        clip_polygon([(0, 0), (1, 0), (1, 1)], [(0, 0), (1, 1)])
