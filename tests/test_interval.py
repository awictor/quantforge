"""Interval arithmetic: rigorous enclosures."""

import math
import random

import pytest

from quantforge import Interval


def _sample(iv, rng):
    return rng.uniform(iv.lo, iv.hi)


def test_arithmetic_enclosure():
    rng = random.Random(1)
    for _ in range(5000):
        a = Interval(rng.uniform(-10, 10), rng.uniform(-10, 10))
        b = Interval(rng.uniform(-10, 10), rng.uniform(-10, 10))
        add, sub, mul = a + b, a - b, a * b
        for _ in range(5):
            x, y = _sample(a, rng), _sample(b, rng)
            assert add.contains(x + y)
            assert sub.contains(x - y)
            assert mul.contains(x * y)


def test_division_enclosure():
    rng = random.Random(2)
    for _ in range(3000):
        a = Interval(rng.uniform(-10, 10), rng.uniform(-10, 10))
        lo = rng.uniform(1, 10)
        b = Interval(lo, lo + rng.uniform(0.1, 5))
        d = a / b
        for _ in range(5):
            x, y = _sample(a, rng), _sample(b, rng)
            assert d.contains(x / y)


def test_known_values():
    assert Interval(1, 2) + Interval(3, 4) == Interval(4, 6)
    assert Interval(1, 2) * Interval(-1, 3) == Interval(-2, 6)
    assert Interval(-2, 3) ** 2 == Interval(0, 9)      # straddles zero
    assert Interval(2, 3) ** 2 == Interval(4, 9)
    assert Interval(-3, -2) ** 3 == Interval(-27, -8)


def test_helpers():
    i = Interval(2, 8)
    assert i.width() == 6 and i.midpoint() == 5
    assert i.contains(5) and not i.contains(9)
    assert Interval(2, 8).intersect(Interval(5, 10)) == Interval(5, 8)
    assert Interval(0, 1).intersect(Interval(2, 3)) is None


def test_monotone_functions_enclose():
    rng = random.Random(3)
    for _ in range(2000):
        lo = rng.uniform(0.1, 5)
        iv = Interval(lo, lo + rng.uniform(0, 3))
        x = _sample(iv, rng)
        assert iv.exp().contains(math.exp(x))
        assert iv.log().contains(math.log(x))
        assert iv.sqrt().contains(math.sqrt(x))


def test_validation():
    with pytest.raises(ValueError):
        Interval(1, 2) / Interval(-1, 1)
    with pytest.raises(ValueError):
        Interval(-1, 2).log()
    with pytest.raises(ValueError):
        Interval(1, 2) ** -1
