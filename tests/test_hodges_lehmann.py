"""Hodges-Lehmann robust location and shift."""

import random
import statistics

import pytest

from quantforge import hodges_lehmann_location, hodges_lehmann_shift


def _hl_ref(x):
    n = len(x)
    w = sorted((x[i] + x[j]) / 2 for i in range(n) for j in range(i, n))
    m = len(w)
    return w[m // 2] if m % 2 else (w[m // 2 - 1] + w[m // 2]) / 2


def _shift_ref(x, y):
    d = sorted(yj - xi for xi in x for yj in y)
    m = len(d)
    return d[m // 2] if m % 2 else (d[m // 2 - 1] + d[m // 2]) / 2


def test_location_of_symmetric():
    assert hodges_lehmann_location([1, 2, 3, 4, 5]) == 3.0


def test_location_matches_reference():
    rng = random.Random(1)
    for _ in range(300):
        n = rng.randint(1, 20)
        xs = [rng.gauss(0, 3) for _ in range(n)]
        assert abs(hodges_lehmann_location(xs) - _hl_ref(xs)) < 1e-9


def test_location_robust_to_outlier():
    x = [10, 11, 9, 12, 10, 11, 9, 10, 13, 8]
    xo = x + [1000.0]
    assert statistics.mean(xo) > 90         # mean dragged away
    assert abs(hodges_lehmann_location(xo) - 10.5) < 1.0   # HL barely moves


def test_shift_exact():
    x = [1, 2, 3, 4, 5]
    y = [v + 5 for v in x]
    assert hodges_lehmann_shift(x, y) == 5


def test_shift_matches_reference():
    rng = random.Random(7)
    for _ in range(200):
        n = rng.randint(1, 12)
        m = rng.randint(1, 12)
        xs = [rng.gauss(0, 1) for _ in range(n)]
        ys = [rng.gauss(0.5, 1) for _ in range(m)]
        assert abs(hodges_lehmann_shift(xs, ys) - _shift_ref(xs, ys)) < 1e-9


def test_shift_antisymmetry():
    assert abs(hodges_lehmann_shift([1, 2, 3], [4, 5, 6])
               + hodges_lehmann_shift([4, 5, 6], [1, 2, 3])) < 1e-12


def test_shift_zero_for_same_distribution():
    rng = random.Random(9)
    a = [rng.gauss(0, 1) for _ in range(200)]
    b = [rng.gauss(0, 1) for _ in range(200)]
    assert abs(hodges_lehmann_shift(a, b)) < 0.3


def test_validation():
    with pytest.raises(ValueError):
        hodges_lehmann_location([])
    with pytest.raises(ValueError):
        hodges_lehmann_shift([1], [])
