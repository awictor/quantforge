"""Tests for Catmull-Rom splines: interpolation, continuity, and passing through points."""

import random

import pytest

from quantforge.catmull_rom import catmull_rom_point, catmull_rom_curve


def _close(a, b, tol=1e-9):
    if isinstance(a, (int, float)):
        return abs(a - b) <= tol * max(1, abs(a), abs(b))
    return all(abs(x - y) <= tol * max(1, abs(x), abs(y)) for x, y in zip(a, b))


def test_fuzz_interpolation_and_continuity():
    rng = random.Random(361)
    for _ in range(4000):
        n = rng.randint(2, 8)
        dim = rng.choice([0, 2, 3])
        if dim == 0:
            pts = [rng.uniform(-10, 10) for _ in range(n)]
        else:
            pts = [tuple(rng.uniform(-10, 10) for _ in range(dim)) for _ in range(n)]
        alpha = rng.choice([0.0, 0.5, 1.0])
        for seg in range(n - 1):
            assert _close(catmull_rom_point(pts, seg, 0.0, alpha), pts[seg])
            assert _close(catmull_rom_point(pts, seg, 1.0, alpha), pts[seg + 1])
        for seg in range(n - 2):
            end = catmull_rom_point(pts, seg, 1.0, alpha)
            start = catmull_rom_point(pts, seg + 1, 0.0, alpha)
            assert _close(end, start)


def test_curve_passes_through_all_points():
    rng = random.Random(362)
    for _ in range(1000):
        n = rng.randint(2, 6)
        pts = [tuple(rng.uniform(-5, 5) for _ in range(2)) for _ in range(n)]
        curve = catmull_rom_curve(pts, samples_per_segment=8)
        for p in pts:
            assert any(_close(p, c) for c in curve)
        assert _close(curve[0], pts[0])
        assert _close(curve[-1], pts[-1])


def test_interpolates_control_points_2d():
    pts = [(0, 0), (1, 1), (2, 0), (3, 1)]
    assert _close(catmull_rom_point(pts, 1, 0.0), (1, 1))
    assert _close(catmull_rom_point(pts, 1, 1.0), (2, 0))


def test_two_point_segment_midpoint():
    assert _close(catmull_rom_point([0.0, 10.0], 0, 0.5), 5.0)


def test_alpha_variants_all_interpolate():
    pts = [(0, 0), (1, 3), (4, 2), (5, 5)]
    for alpha in (0.0, 0.5, 1.0):
        for seg in range(3):
            assert _close(catmull_rom_point(pts, seg, 0.0, alpha), pts[seg])
            assert _close(catmull_rom_point(pts, seg, 1.0, alpha), pts[seg + 1])


def test_coincident_points_do_not_crash():
    pts = [(0, 0), (0, 0), (1, 1)]
    v = catmull_rom_point(pts, 1, 0.5)
    assert v is not None


def test_curve_sample_count():
    pts = [0, 1, 2, 3]  # 3 segments
    curve = catmull_rom_curve(pts, samples_per_segment=5)
    # first point + (samples-1) per segment
    assert len(curve) == 1 + 3 * (5 - 1)


def test_straight_line_stays_straight():
    # collinear points -> the spline should stay near the line at the midpoint
    pts = [(0, 0), (1, 1), (2, 2), (3, 3)]
    mid = catmull_rom_point(pts, 1, 0.5)
    assert _close(mid, (1.5, 1.5), tol=1e-6)


def test_too_few_points_raises():
    with pytest.raises(ValueError):
        catmull_rom_point([1], 0, 0.5)
    with pytest.raises(ValueError):
        catmull_rom_curve([1], 4)


def test_bad_segment_raises():
    with pytest.raises(IndexError):
        catmull_rom_point([1, 2], 5, 0.5)


def test_too_few_samples_raises():
    with pytest.raises(ValueError):
        catmull_rom_curve([1, 2], 1)
