"""Partial and semi-partial correlation."""

import math
import random

import pytest

from quantforge import partial_correlation, semipartial_correlation
from quantforge.partial_corr import _pearson


def test_spurious_correlation_drops_out():
    rng = random.Random(1)
    n = 2000
    z = [rng.gauss(0, 1) for _ in range(n)]
    x = [z[i] + rng.gauss(0, 0.3) for i in range(n)]
    y = [z[i] + rng.gauss(0, 0.3) for i in range(n)]
    assert _pearson(x, y) > 0.8              # strong raw correlation
    assert abs(partial_correlation(x, y, z)) < 0.1   # vanishes controlling Z


def test_matches_three_variable_formula():
    rng = random.Random(2)
    for _ in range(300):
        z = [rng.gauss(0, 1) for _ in range(200)]
        x = [z[i] * 0.5 + rng.gauss(0, 1) for i in range(200)]
        y = [z[i] * 0.3 + x[i] * 0.2 + rng.gauss(0, 1) for i in range(200)]
        rxy, rxz, ryz = _pearson(x, y), _pearson(x, z), _pearson(y, z)
        formula = (rxy - rxz * ryz) / math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
        assert abs(partial_correlation(x, y, z) - formula) < 1e-9


def test_uncorrelated_control_preserves_correlation():
    rng = random.Random(3)
    z = [rng.gauss(0, 1) for _ in range(2000)]
    x = [rng.gauss(0, 1) for _ in range(2000)]
    y = [x[i] * 0.7 + rng.gauss(0, 0.5) for i in range(2000)]
    assert abs(partial_correlation(x, y, z) - _pearson(x, y)) < 0.02


def test_multiple_controls():
    rng = random.Random(4)
    z1 = [rng.gauss(0, 1) for _ in range(1000)]
    z2 = [rng.gauss(0, 1) for _ in range(1000)]
    x = [z1[i] + z2[i] + rng.gauss(0, 0.3) for i in range(1000)]
    y = [z1[i] + z2[i] + rng.gauss(0, 0.3) for i in range(1000)]
    assert abs(partial_correlation(x, y, [z1, z2])) < 0.1


def test_semipartial_bounded():
    rng = random.Random(5)
    z = [rng.gauss(0, 1) for _ in range(500)]
    x = [rng.gauss(0, 1) for _ in range(500)]
    y = [z[i] + x[i] * 0.3 + rng.gauss(0, 0.5) for i in range(500)]
    sp = semipartial_correlation(x, y, z)
    assert -1.0 <= sp <= 1.0


def test_validation():
    with pytest.raises(ValueError):
        partial_correlation([1, 1, 1], [1, 2, 3], [1, 2, 3])
