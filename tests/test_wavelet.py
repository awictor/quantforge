"""Haar discrete wavelet transform."""

import math
import random

import pytest

from quantforge import haar_dwt, haar_idwt, wavelet_energy


def test_perfect_reconstruction():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(16)]
    a, d = haar_dwt(x)
    xr = haar_idwt(a, d)
    assert max(abs(xr[i] - x[i]) for i in range(16)) < 1e-10


def test_constant_has_no_detail():
    a, d = haar_dwt([5.0] * 8)
    assert all(abs(v) < 1e-12 for lvl in d for v in lvl)


def test_parseval_energy_preserved():
    rng = random.Random(2)
    x = [rng.gauss(0, 1) for _ in range(16)]
    a, d = haar_dwt(x)
    coef_energy = sum(v * v for v in a) + sum(v * v for lvl in d for v in lvl)
    assert abs(coef_energy - sum(v * v for v in x)) < 1e-9


def test_energy_fractions_sum_to_one():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(32)]
    e = wavelet_energy(x)
    assert abs(sum(e["detail"]) + e["approx"] - 1.0) < 1e-9


def test_constant_energy_in_approximation():
    e = wavelet_energy([3.0] * 8)
    assert abs(e["approx"] - 1.0) < 1e-9


def test_single_level_average_difference():
    a, d = haar_dwt([3.0, 1.0], levels=1)
    assert abs(a[0] - 4.0 / math.sqrt(2)) < 1e-12
    assert abs(d[0][0] - 2.0 / math.sqrt(2)) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        haar_dwt([1.0, 2.0, 3.0])          # not a power of two
    with pytest.raises(ValueError):
        haar_dwt([1.0, 2.0], levels=5)     # too many levels
    with pytest.raises(ValueError):
        wavelet_energy([0.0, 0.0])          # zero energy
