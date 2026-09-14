"""Modular arithmetic: extended GCD, inverse, CRT, discrete log."""

import math
import random

import pytest

from quantforge import (
    extended_gcd,
    mod_inverse,
    chinese_remainder,
    mod_pow,
    discrete_log,
)


def test_extended_gcd_bezout():
    rng = random.Random(1)
    for _ in range(1000):
        a = rng.randint(1, 10 ** 6)
        b = rng.randint(1, 10 ** 6)
        g, x, y = extended_gcd(a, b)
        assert a * x + b * y == g
        assert g == math.gcd(a, b)


def test_mod_inverse():
    rng = random.Random(2)
    for _ in range(1000):
        m = rng.randint(2, 10 ** 6)
        a = rng.randint(1, m - 1)
        if math.gcd(a, m) == 1:
            inv = mod_inverse(a, m)
            assert (a * inv) % m == 1
            assert 0 <= inv < m


def test_mod_inverse_requires_coprime():
    with pytest.raises(ValueError):
        mod_inverse(4, 8)


def test_crt_reconstructs():
    x_true = 1234567
    moduli = [7, 11, 13, 17, 19]
    rem = [x_true % m for m in moduli]
    x, M = chinese_remainder(rem, moduli)
    assert M == 7 * 11 * 13 * 17 * 19
    assert x == x_true % M
    assert chinese_remainder([2, 3, 2], [3, 5, 7]) == (23, 105)


def test_crt_requires_coprime():
    with pytest.raises(ValueError):
        chinese_remainder([1, 2], [4, 6])


def test_mod_pow_including_negative():
    assert mod_pow(3, 4, 5) == 1
    assert mod_pow(3, -1, 5) == 2          # 3*2 = 6 = 1 (mod 5)
    assert mod_pow(2, 10, 1000) == 24


def test_discrete_log_prime_field():
    p, g = 1019, 2
    rng = random.Random(3)
    for _ in range(50):
        x = rng.randint(0, p - 2)
        t = pow(g, x, p)
        xr = discrete_log(g, t, p)
        assert xr is not None and pow(g, xr, p) == t
    assert discrete_log(3, 13, 17) == 4
    assert discrete_log(2, 0, 7) is None    # no solution


def test_validation():
    with pytest.raises(ValueError):
        mod_inverse(3, 0)
    with pytest.raises(ValueError):
        discrete_log(2, 3, 0)
    with pytest.raises(ValueError):
        chinese_remainder([1], [])
