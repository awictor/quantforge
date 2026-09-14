"""Tests for prime-field polynomial arithmetic and Lagrange interpolation."""

import random

import pytest

from quantforge.poly_mod import (
    poly_eval_mod,
    poly_add_mod,
    poly_mul_mod,
    lagrange_interpolate_mod,
)


def _brute_mul(a, b, mod):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = (out[i + j] + ai * bj) % mod
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def test_fuzz_interpolation_and_arithmetic():
    rng = random.Random(281)
    for mod in (997, 998244353):
        for _ in range(2000):
            deg = rng.randint(0, 6)
            coeffs = [rng.randint(0, mod - 1) for _ in range(deg + 1)]
            while len(coeffs) > 1 and coeffs[-1] == 0:
                coeffs.pop()
            n = len(coeffs)
            xs = rng.sample(range(mod), n)
            pts = [(x, poly_eval_mod(coeffs, x, mod)) for x in xs]
            rec = lagrange_interpolate_mod(pts, mod)
            for x, y in pts:
                assert poly_eval_mod(rec, x, mod) == y
            assert rec == coeffs
            a = [rng.randint(0, mod - 1) for _ in range(rng.randint(1, 5))]
            b = [rng.randint(0, mod - 1) for _ in range(rng.randint(1, 5))]
            assert poly_mul_mod(a, b, mod) == _brute_mul(a, b, mod)
            x = rng.randint(0, mod - 1)
            assert poly_eval_mod(poly_add_mod(a, b, mod), x, mod) == (
                poly_eval_mod(a, x, mod) + poly_eval_mod(b, x, mod)
            ) % mod


def test_interpolate_quadratic():
    # f(0)=1, f(1)=3, f(2)=7 -> x^2 + x + 1
    assert lagrange_interpolate_mod([(0, 1), (1, 3), (2, 7)], 997) == [1, 1, 1]


def test_interpolate_single_point():
    assert lagrange_interpolate_mod([(5, 42)], 997) == [42]


def test_interpolate_line():
    # f(0)=2, f(1)=5 -> 3x + 2
    assert lagrange_interpolate_mod([(0, 2), (1, 5)], 997) == [2, 3]


def test_poly_eval_horner():
    assert poly_eval_mod([1, 2, 3], 2, 997) == 17  # 1 + 2*2 + 3*4
    assert poly_eval_mod([], 5, 997) == 0
    assert poly_eval_mod([7], 100, 997) == 7


def test_poly_add():
    assert poly_add_mod([1, 2, 3], [4, 5], 997) == [5, 7, 3]
    assert poly_add_mod([1], [996], 997) == [0]  # wraps to 0 mod 997


def test_poly_mul():
    # (1 + x)(1 + x) = 1 + 2x + x^2
    assert poly_mul_mod([1, 1], [1, 1], 997) == [1, 2, 1]
    assert poly_mul_mod([], [1, 2], 997) == []


def test_multiplication_wraps_mod():
    assert poly_mul_mod([500], [500], 997) == [500 * 500 % 997]


def test_interpolation_recovers_zero_polynomial():
    pts = [(0, 0), (1, 0), (2, 0)]
    assert lagrange_interpolate_mod(pts, 997) == [0]


def test_duplicate_x_raises():
    with pytest.raises(ValueError):
        lagrange_interpolate_mod([(1, 2), (1, 3)], 997)


def test_duplicate_x_modular_raises():
    # 1 and 1+mod collapse to the same residue
    with pytest.raises(ValueError):
        lagrange_interpolate_mod([(1, 2), (1 + 997, 3)], 997)


def test_shamir_secret_sharing_roundtrip():
    # secret = constant term; interpolate from any t shares recovers it
    mod = 998244353
    secret = 123456
    coeffs = [secret, 111, 222]  # degree-2 -> needs 3 shares
    shares = [(x, poly_eval_mod(coeffs, x, mod)) for x in (1, 2, 3, 4, 5)]
    rec = lagrange_interpolate_mod(shares[:3], mod)
    assert rec[0] == secret
    assert lagrange_interpolate_mod(shares[2:5], mod)[0] == secret
