"""Directional (circular) statistics for angles."""

import math
import random

import pytest

from quantforge import (
    circular_mean,
    resultant_length,
    circular_variance,
    circular_std,
    rayleigh_test,
)

d2r = math.radians


def test_wraparound_mean():
    assert abs(math.degrees(circular_mean([d2r(350), d2r(10)]))) < 1e-9    # -> 0
    assert abs(math.degrees(circular_mean([d2r(10), d2r(20), d2r(30)])) - 20) < 1e-9


def test_resultant_length_extremes():
    assert abs(resultant_length([1.0, 1.0, 1.0]) - 1.0) < 1e-12
    uniform = [2 * math.pi * k / 360 for k in range(360)]
    assert resultant_length(uniform) < 1e-6
    assert circular_variance([1.0, 1.0]) == 0.0
    assert circular_std([0.5, 0.5]) == 0.0


def test_cancelling_vectors_raise():
    with pytest.raises(ValueError):
        circular_mean([0.0, math.pi])


def test_rayleigh_rejects_concentrated():
    rng = random.Random(1)
    conc = [d2r(90) + rng.gauss(0, 0.2) for _ in range(50)]
    R, p = rayleigh_test(conc)
    assert R > 0.9 and p < 0.01
    assert abs(math.degrees(circular_mean(conc)) - 90) < 5


def test_rayleigh_accepts_uniform():
    rng = random.Random(2)
    unif = [rng.uniform(0, 2 * math.pi) for _ in range(200)]
    R, p = rayleigh_test(unif)
    assert R < 0.15 and p > 0.05


def test_validation():
    with pytest.raises(ValueError):
        resultant_length([])
    with pytest.raises(ValueError):
        rayleigh_test([1.0])
