"""Gaver-Stehfest numerical Laplace inversion."""

import math

import pytest

from quantforge import laplace_inversion


def test_known_pairs_smooth():
    cases = [
        (lambda s: 1.0 / s, lambda t: 1.0),
        (lambda s: 1.0 / (s * s), lambda t: t),
        (lambda s: 2.0 / s ** 3, lambda t: t * t),
        (lambda s: 1.0 / (s + 2.0), lambda t: math.exp(-2 * t)),
        (lambda s: math.sqrt(math.pi) / (2 * s ** 1.5), lambda t: math.sqrt(t)),
    ]
    for F, f in cases:
        for t in (0.5, 1.0, 2.0):
            got = laplace_inversion(F, t, N=12)
            assert abs(got - f(t)) < 0.01 * (abs(f(t)) + 0.01)


def test_decaying_exponential():
    F = lambda s: 1.0 / (s + 2.0)
    assert abs(laplace_inversion(F, 1.0) - math.exp(-2.0)) < 1e-3


def test_cosine_short_time():
    F = lambda s: s / (s * s + 1.0)
    assert abs(laplace_inversion(F, 0.5) - math.cos(0.5)) < 1e-3


def test_ramp_scales_linearly():
    F = lambda s: 1.0 / (s * s)
    assert abs(laplace_inversion(F, 3.0) - 3.0) < 1e-3


def test_validation():
    with pytest.raises(ValueError):
        laplace_inversion(lambda s: 1.0 / s, 0.0)
    with pytest.raises(ValueError):
        laplace_inversion(lambda s: 1.0 / s, 1.0, N=11)
