"""One-dimensional minimizers: golden-section and Brent."""

import math

import pytest

from quantforge import golden_section_min, brent_min


def test_quadratic_minimum():
    f = lambda x: (x - 3) ** 2 + 1
    xg, fg = golden_section_min(f, -10, 10)
    xb, fb = brent_min(f, -10, 10)
    assert abs(xg - 3) < 1e-6 and abs(fg - 1) < 1e-6
    assert abs(xb - 3) < 1e-8 and abs(fb - 1) < 1e-8


def test_cosine_minimum_at_pi():
    x, _ = brent_min(math.cos, 0, 2 * math.pi)
    assert abs(x - math.pi) < 1e-6


def test_quartic_minimum_at_zero():
    x, _ = brent_min(lambda x: x ** 4, -2, 2)
    assert abs(x) < 1e-3


def test_gaussian_peak_via_golden():
    x, _ = golden_section_min(lambda x: -math.exp(-(x - 1) ** 2), -5, 5)
    assert abs(x - 1) < 1e-6


def test_methods_agree():
    f = lambda x: (x - 3) ** 2 + 1
    xg, _ = golden_section_min(f, -10, 10)
    xb, _ = brent_min(f, -10, 10)
    assert abs(xg - xb) < 1e-5


def test_validation():
    f = lambda x: x * x
    with pytest.raises(ValueError):
        golden_section_min(f, 5, 1)
    with pytest.raises(ValueError):
        brent_min(f, 5, 1)
