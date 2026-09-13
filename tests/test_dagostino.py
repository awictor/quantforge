"""D'Agostino-Pearson K^2 omnibus normality test."""

import random

import pytest

from quantforge import dagostino_k2


def test_null_calibration():
    rng = random.Random(3)
    rej = 0
    N = 1000
    for _ in range(N):
        d = [rng.gauss(0, 1) for _ in range(100)]
        if dagostino_k2(d)["p_value"] < 0.05:
            rej += 1
    assert 0.03 <= rej / N <= 0.07


def test_rejects_skew():
    rng = random.Random(5)
    d = [rng.expovariate(1.0) for _ in range(200)]
    r = dagostino_k2(d)
    assert r["p_value"] < 0.01
    assert r["z_skew"] > 3


def test_rejects_heavy_tails():
    rng = random.Random(7)
    d = [rng.gauss(0, 1) * (3 if rng.random() < 0.1 else 1) for _ in range(500)]
    r = dagostino_k2(d)
    assert r["p_value"] < 0.01
    assert r["z_kurt"] > 2


def test_light_tails_negative_kurtosis_z():
    rng = random.Random(9)
    d = [rng.uniform(-1, 1) for _ in range(500)]
    assert dagostino_k2(d)["z_kurt"] < 0


def test_clean_normal_high_p():
    rng = random.Random(11)
    d = [rng.gauss(5, 2) for _ in range(2000)]
    assert dagostino_k2(d)["p_value"] > 0.05


def test_power_on_exponential():
    rng = random.Random(13)
    detected = sum(1 for _ in range(200)
                   if dagostino_k2([rng.expovariate(1) for _ in range(50)])["p_value"] < 0.05)
    assert detected / 200 > 0.8


def test_validation():
    with pytest.raises(ValueError):
        dagostino_k2([1, 2, 3])
    with pytest.raises(ValueError):
        dagostino_k2([5.0] * 30)             # zero variance
