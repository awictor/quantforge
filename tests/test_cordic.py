"""Tests for CORDIC sin/cos/atan2/hypot, cross-checked against the math module."""

import math
import random

from quantforge.cordic import cordic_sincos, cordic_atan2, cordic_hypot


def _close(a, b, tol=1e-9):
    return abs(a - b) <= tol


def test_fuzz_sincos_full_range():
    rng = random.Random(531)
    for _ in range(10000):
        th = rng.uniform(-4 * math.pi, 4 * math.pi)
        c, s = cordic_sincos(th)
        assert _close(c, math.cos(th), 1e-9)
        assert _close(s, math.sin(th), 1e-9)


def test_fuzz_atan2_and_hypot():
    rng = random.Random(532)
    for _ in range(10000):
        x = rng.uniform(-100, 100)
        y = rng.uniform(-100, 100)
        assert _close(cordic_atan2(y, x), math.atan2(y, x), 1e-9)
        h = math.hypot(x, y)
        assert _close(cordic_hypot(x, y), h, 1e-6 * max(1, abs(h)))


def test_sincos_landmarks():
    for th, (ec, es) in [
        (0.0, (1, 0)),
        (math.pi / 2, (0, 1)),
        (math.pi, (-1, 0)),
        (-math.pi / 2, (0, -1)),
        (math.pi / 4, (math.sqrt(2) / 2, math.sqrt(2) / 2)),
    ]:
        c, s = cordic_sincos(th)
        assert _close(c, ec, 1e-9)
        assert _close(s, es, 1e-9)


def test_atan2_quadrants():
    assert _close(cordic_atan2(1, 1), math.pi / 4)
    assert _close(cordic_atan2(1, -1), 3 * math.pi / 4)
    assert _close(cordic_atan2(-1, -1), -3 * math.pi / 4)
    assert _close(cordic_atan2(-1, 1), -math.pi / 4)
    assert _close(cordic_atan2(1, 0), math.pi / 2)
    assert _close(cordic_atan2(-1, 0), -math.pi / 2)
    assert _close(cordic_atan2(0, -1), math.pi)


def test_hypot_known():
    assert _close(cordic_hypot(3, 4), 5.0, 1e-6)
    assert _close(cordic_hypot(5, 12), 13.0, 1e-6)
    assert _close(cordic_hypot(-3, -4), 5.0, 1e-6)


def test_zero_cases():
    assert cordic_atan2(0, 0) == 0.0
    assert cordic_hypot(0, 0) == 0.0
    c, s = cordic_sincos(0.0)
    assert _close(c, 1.0) and _close(s, 0.0)


def test_periodicity():
    c1, s1 = cordic_sincos(0.5)
    c2, s2 = cordic_sincos(0.5 + 2 * math.pi)
    assert _close(c1, c2, 1e-9)
    assert _close(s1, s2, 1e-9)


def test_pythagorean_identity():
    rng = random.Random(533)
    for _ in range(1000):
        th = rng.uniform(-10, 10)
        c, s = cordic_sincos(th)
        assert _close(c * c + s * s, 1.0, 1e-9)


def test_axis_hypot():
    assert _close(cordic_hypot(7, 0), 7.0, 1e-6)
    assert _close(cordic_hypot(0, 9), 9.0, 1e-6)
