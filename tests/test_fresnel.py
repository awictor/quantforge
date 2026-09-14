"""Tests for Fresnel integrals and the Dawson function, vs numerical integration."""

import math
import random

import pytest

from quantforge.fresnel import fresnel_c, fresnel_s, dawson


def _close(a, b, tol=1e-5):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _num_c(x, N=200000):
    if x == 0:
        return 0.0
    h = x / N
    s = 0.0
    for i in range(N + 1):
        t = i * h
        w = 0.5 if i in (0, N) else 1.0
        s += w * math.cos(math.pi / 2 * t * t)
    return s * h


def _num_s(x, N=200000):
    if x == 0:
        return 0.0
    h = x / N
    s = 0.0
    for i in range(N + 1):
        t = i * h
        w = 0.5 if i in (0, N) else 1.0
        s += w * math.sin(math.pi / 2 * t * t)
    return s * h


def _num_d(x, N=200000):
    if x == 0:
        return 0.0
    h = x / N
    s = 0.0
    for i in range(N + 1):
        t = i * h
        w = 0.5 if i in (0, N) else 1.0
        s += w * math.exp(t * t)
    return math.exp(-x * x) * s * h


def test_zero():
    assert _close(fresnel_c(0), 0.0, 1e-9)
    assert _close(fresnel_s(0), 0.0, 1e-9)
    assert _close(dawson(0), 0.0, 1e-9)


def test_fresnel_vs_integration_fast():
    for x in (0.3, 1.0, 2.5, 3.5):
        assert _close(fresnel_c(x), _num_c(x, 40000), 1e-4)
        assert _close(fresnel_s(x), _num_s(x, 40000), 1e-4)


@pytest.mark.slow
def test_fuzz_fresnel_vs_integration():
    rng = random.Random(571)
    for _ in range(300):
        x = rng.uniform(0.05, 7)
        assert _close(fresnel_c(x), _num_c(x), 1e-5)
        assert _close(fresnel_s(x), _num_s(x), 1e-5)


def test_fresnel_odd_symmetry():
    for x in (0.5, 1.3, 3.0, 5.0):
        assert _close(fresnel_c(-x), -fresnel_c(x))
        assert _close(fresnel_s(-x), -fresnel_s(x))


@pytest.mark.slow
def test_fresnel_asymptotic_region():
    for x in (4.0, 5.0, 8.0):
        assert _close(fresnel_c(x), _num_c(x, 800000), 1e-5)
        assert _close(fresnel_s(x), _num_s(x, 800000), 1e-5)


def test_fresnel_approaches_half():
    # oscillates around 0.5 with amplitude ~ 1/(pi x), so use a large x and loose tol
    assert abs(fresnel_c(2000) - 0.5) < 1e-3
    assert abs(fresnel_s(2000) - 0.5) < 1e-3


def test_dawson_known_peak():
    # D peaks at x ~ 0.9241 with value ~ 0.5410
    assert _close(dawson(0.9241388730), 0.5410442246, 1e-6)


def test_dawson_vs_integration():
    for x in (0.3, 0.5, 1.0, 2.0, 4.0):
        assert _close(dawson(x), _num_d(x), 1e-5)


def test_dawson_odd_symmetry():
    for x in (0.5, 1.5, 3.0):
        assert _close(dawson(-x), -dawson(x))


def test_dawson_large_x_decay():
    # D(x) ~ 1/(2x) for large x
    assert _close(dawson(10.0), 1 / 20.0, 1e-3)
    assert _close(dawson(50.0), 1 / 100.0, 1e-4)


def test_fresnel_small_x_matches_leading_term():
    # C(x) ~ x, S(x) ~ pi x^3 / 6 for small x
    x = 0.01
    assert _close(fresnel_c(x), x, 1e-6)
    assert _close(fresnel_s(x), math.pi * x ** 3 / 6, 1e-6)
