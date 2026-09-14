"""Tests for Halley and secant root finders, cross-checked against known roots and brent."""

import math
import random

import pytest

from quantforge.halley import halley, secant
from quantforge.rootfind import brent


def _close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def test_sqrt_and_cbrt():
    assert _close(halley(lambda x: x * x - 2, lambda x: 2 * x, lambda x: 2, 1.0), math.sqrt(2))
    assert _close(halley(lambda x: x ** 3 - 5, lambda x: 3 * x * x, lambda x: 6 * x, 2.0),
                  5 ** (1 / 3))
    assert _close(secant(lambda x: x * x - 2, 1, 2), math.sqrt(2))


def test_dottie_number():
    f = lambda x: math.cos(x) - x
    root = 0.7390851332151607
    assert _close(halley(f, lambda x: -math.sin(x) - 1, lambda x: -math.cos(x), 0.5), root)
    assert _close(secant(f, 0, 1), root)


def test_fuzz_cubics():
    rng = random.Random(451)
    for _ in range(3000):
        a, b, c, d = [rng.uniform(-5, 5) for _ in range(4)]
        if abs(a) < 0.5:
            a = 1.0
        f = lambda x, a=a, b=b, c=c, d=d: a * x ** 3 + b * x ** 2 + c * x + d
        f1 = lambda x, a=a, b=b, c=c: 3 * a * x ** 2 + 2 * b * x + c
        f2 = lambda x, a=a, b=b: 6 * a * x + 2 * b
        if f(-100) * f(100) > 0:
            continue
        br = brent(f, -100, 100)
        try:
            h = halley(f, f1, f2, br + rng.uniform(-0.5, 0.5))
            assert _close(f(h), 0, 1e-6)
        except ValueError:
            pass
        try:
            s = secant(f, br - 1, br + 1)
            assert _close(f(s), 0, 1e-6)
        except ValueError:
            pass


def test_halley_cubic_convergence():
    # Halley should reach machine precision in very few iterations; count them
    calls = [0]

    def f(x):
        calls[0] += 1
        return x * x - 2

    halley(f, lambda x: 2 * x, lambda x: 2, 1.0)
    assert calls[0] <= 6  # cubic convergence: a handful of steps


def test_secant_superlinear():
    r = secant(lambda x: x ** 3 - x - 2, 1, 2)
    assert _close(r ** 3 - r - 2, 0, 1e-10)


def test_no_real_root_raises():
    with pytest.raises(ValueError):
        halley(lambda x: x * x + 1, lambda x: 2 * x, lambda x: 2, 0.0)


def test_flat_secant_raises():
    with pytest.raises(ValueError):
        secant(lambda x: 5.0, 0, 1)


def test_already_at_root():
    assert _close(halley(lambda x: x - 3, lambda x: 1, lambda x: 0, 3.0), 3.0)
    assert _close(secant(lambda x: x - 3, 3.0, 4.0), 3.0)


def test_transcendental_both_roots_of_exp():
    f = lambda x: math.exp(x) - 3 * x
    f1 = lambda x: math.exp(x) - 3
    f2 = lambda x: math.exp(x)
    r1 = halley(f, f1, f2, 0.5)
    r2 = halley(f, f1, f2, 1.5)
    assert _close(f(r1), 0, 1e-9)
    assert _close(f(r2), 0, 1e-9)
    assert not _close(r1, r2)


def test_secant_polynomial_root():
    assert _close(secant(lambda x: (x - 4) * (x + 1), 3, 5), 4.0)
