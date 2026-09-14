"""Interpolation on a 2-D regular grid: bilinear and nearest."""

import random

import pytest

from quantforge import bilinear_interp, nearest_interp

XS = [0, 1, 2, 3]
YS = [0, 1, 2]


def _planar():
    return [[2 * x + 3 * y + 1 for x in XS] for y in YS]


def test_exact_at_nodes():
    z = _planar()
    for j, y in enumerate(YS):
        for i, x in enumerate(XS):
            assert abs(bilinear_interp(XS, YS, z, x, y) - z[j][i]) < 1e-12


def test_planar_field_is_exact():
    z = _planar()
    rng = random.Random(1)
    for _ in range(1000):
        x, y = rng.uniform(0, 3), rng.uniform(0, 2)
        assert abs(bilinear_interp(XS, YS, z, x, y) - (2 * x + 3 * y + 1)) < 1e-9


def test_cell_center_is_corner_mean():
    z = [[1, 2, 4, 8], [3, 5, 9, 17], [7, 11, 19, 35]]
    c = bilinear_interp([0, 1, 2, 3], [0, 1, 2], z, 0.5, 0.5)
    assert abs(c - (1 + 2 + 3 + 5) / 4) < 1e-12


def test_clamp_out_of_range():
    z = _planar()
    assert bilinear_interp(XS, YS, z, -5, -5) == z[0][0]
    assert bilinear_interp(XS, YS, z, 99, 99) == z[-1][-1]


def test_nearest_snaps():
    z = _planar()
    assert nearest_interp(XS, YS, z, 0.4, 0.4) == z[0][0]
    assert nearest_interp(XS, YS, z, 0.6, 0.6) == z[1][1]


def test_validation():
    with pytest.raises(ValueError):
        bilinear_interp([0], [0, 1], [[1], [2]], 0, 0)
    with pytest.raises(ValueError):
        bilinear_interp([0, 1], [0, 1], [[1, 2]], 0, 0)
