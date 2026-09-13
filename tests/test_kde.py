"""Gaussian kernel density estimation."""

import math
import random

import pytest

from quantforge import kde, silverman_bandwidth, scott_bandwidth, kde_function


def test_integrates_to_one():
    rng = random.Random(3)
    data = [rng.gauss(0, 1) for _ in range(500)]
    xs = [i * 0.05 for i in range(-160, 161)]
    vals = kde(data, xs)
    integral = sum((vals[i] + vals[i + 1]) / 2 * 0.05 for i in range(len(xs) - 1))
    assert abs(integral - 1.0) < 0.02


def test_non_negative():
    rng = random.Random(3)
    data = [rng.gauss(0, 1) for _ in range(200)]
    assert all(v >= 0 for v in kde(data, [i * 0.1 for i in range(-50, 50)]))


def test_recovers_normal_density():
    rng = random.Random(5)
    big = [rng.gauss(0, 1) for _ in range(20000)]
    assert abs(kde(big, 0.0) - 1 / math.sqrt(2 * math.pi)) < 0.02


def test_scaled_normal():
    rng = random.Random(7)
    d = [rng.gauss(5, 2) for _ in range(20000)]
    assert abs(kde(d, 5.0) - 1 / (2 * math.sqrt(2 * math.pi))) < 0.01


def test_bimodal_peaks_above_valley():
    rng = random.Random(9)
    bim = ([rng.gauss(-3, 0.5) for _ in range(5000)]
           + [rng.gauss(3, 0.5) for _ in range(5000)])
    f = kde_function(bim)
    assert f(-3) > f(0) and f(3) > f(0)


def test_bandwidths_positive():
    rng = random.Random(1)
    data = [rng.gauss(0, 1) for _ in range(500)]
    assert silverman_bandwidth(data) > 0
    assert scott_bandwidth(data) > 0


def test_scalar_list_consistency():
    rng = random.Random(2)
    data = [rng.gauss(0, 1) for _ in range(100)]
    assert kde(data, 0.5) == kde(data, [0.5])[0]


def test_validation():
    with pytest.raises(ValueError):
        kde([1.0], 0.0)
    with pytest.raises(ValueError):
        kde([1.0, 2.0], 0.0, bandwidth=-1)
    with pytest.raises(ValueError):
        silverman_bandwidth([1.0])
