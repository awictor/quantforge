"""Autocorrelation and partial-autocorrelation functions."""

import random

import pytest

from quantforge import acf, pacf


def _ar1(n, phi, seed):
    rng = random.Random(seed)
    y = [0.0]
    for _ in range(n):
        y.append(phi * y[-1] + rng.gauss(0, 1))
    return y[1:]


def _ma1(n, theta, seed):
    rng = random.Random(seed)
    e = [rng.gauss(0, 1) for _ in range(n + 1)]
    return [e[t] + theta * e[t - 1] for t in range(1, n + 1)]


def test_acf_lag_zero_is_one():
    assert acf(_ar1(1000, 0.6, 1), 10)[0] == 1.0


def test_ar1_acf_decays_geometrically():
    phi = 0.6
    a = acf(_ar1(20000, phi, 1), 5)
    assert abs(a[1] - phi) < 0.03
    assert abs(a[2] - phi ** 2) < 0.04
    assert abs(a[3] - phi ** 3) < 0.05


def test_ar1_pacf_cuts_off_after_lag_one():
    phi = 0.6
    p = pacf(_ar1(20000, phi, 1), 10)
    assert abs(p[1] - phi) < 0.03
    assert abs(p[2]) < 0.05
    assert abs(p[3]) < 0.05


def test_ma1_acf_cuts_off_after_lag_one():
    theta = 0.7
    a = acf(_ma1(20000, theta, 2), 10)
    assert abs(a[1] - theta / (1 + theta ** 2)) < 0.03
    assert abs(a[2]) < 0.05


def test_white_noise_acf_small():
    rng = random.Random(3)
    wn = [rng.gauss(0, 1) for _ in range(5000)]
    a = acf(wn, 10)
    assert max(abs(v) for v in a[1:]) < 0.05


def test_validation():
    wn = [0.1, 0.2, 0.3, 0.4]
    with pytest.raises(ValueError):
        acf([1.0])
    with pytest.raises(ValueError):
        acf(wn, 10)               # nlags >= n
    with pytest.raises(ValueError):
        pacf(wn, 0)               # nlags < 1
