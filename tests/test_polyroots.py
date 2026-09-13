"""All-roots polynomial solver (Durand-Kerner)."""

import cmath
import random

import pytest

from quantforge import polynomial_roots


def _all_matched(roots, expected, tol=1e-6):
    return all(min(abs(r - e) for r in roots) < tol for e in expected)


def test_quadratic():
    roots = polynomial_roots([1, -3, 2])       # (x-1)(x-2)
    assert _all_matched(roots, [1, 2])


def test_complex_roots():
    roots = polynomial_roots([1, 0, 1])        # x^2 + 1
    assert _all_matched(roots, [1j, -1j])


def test_cubic():
    roots = polynomial_roots([1, -6, 11, -6])  # (x-1)(x-2)(x-3)
    assert _all_matched(roots, [1, 2, 3])


def test_roots_of_unity():
    roots = polynomial_roots([1, 0, 0, 0, 0, -1])   # x^5 - 1
    expected = [cmath.exp(2j * cmath.pi * k / 5) for k in range(5)]
    assert _all_matched(roots, expected)


def test_double_root():
    roots = polynomial_roots([1, -2, 1])       # (x-1)^2
    assert all(abs(r - 1) < 1e-4 for r in roots)


def test_random_polynomials():
    rng = random.Random(3)
    for _ in range(200):
        n = rng.randint(1, 6)
        rts = [complex(rng.uniform(-5, 5), rng.uniform(-5, 5)) for _ in range(n)]
        coeffs = [1 + 0j]
        for r in rts:
            nc = [0j] * (len(coeffs) + 1)
            for i, cc in enumerate(coeffs):
                nc[i] += cc
                nc[i + 1] -= cc * r
            coeffs = nc
        assert _all_matched(polynomial_roots(coeffs), rts, tol=1e-5)


def test_leading_coefficient_and_residuals():
    coeffs = [2, -3, 0, 1, -5]
    mc = [c / 2 for c in coeffs]
    for r in polynomial_roots(coeffs):
        val = 0j
        for c in mc:
            val = val * r + c
        assert abs(val) < 1e-8


def test_validation():
    with pytest.raises(ValueError):
        polynomial_roots([5])                  # degree 0
    with pytest.raises(ValueError):
        polynomial_roots([0, 0])               # all zero
