"""Dense polynomial arithmetic."""

import random

import pytest

from quantforge import (poly_add, poly_sub, poly_mul, poly_divmod,
                        poly_derivative, poly_integral, poly_eval, poly_gcd)


def test_add_sub():
    assert poly_add([1, 2, 3], [4, 5]) == [5.0, 7.0, 3.0]
    assert poly_sub([1, 2, 3], [1, 2, 3]) == [0.0]


def test_mul_matches_convolution():
    from quantforge import convolve
    a, b = [1, 2, 3], [4, 5, 6, 7]
    pm = poly_mul(a, b)
    cv = list(convolve(a, b))          # FFT-based, carries tiny float noise
    assert len(pm) == len(cv)
    assert all(abs(pm[i] - cv[i]) < 1e-9 for i in range(len(pm)))


def test_divmod_exact():
    q, r = poly_divmod([-1, 0, 1], [-1, 1])       # (x^2-1)/(x-1)
    assert q == [1.0, 1.0]
    assert r == [0.0]


def test_divmod_identity_via_evaluation():
    rng = random.Random(9)
    for _ in range(500):
        n = [rng.uniform(-5, 5) for _ in range(rng.randint(1, 8))]
        d = [rng.uniform(-5, 5) for _ in range(rng.randint(1, 5))]
        if d[-1] == 0:
            continue
        q, r = poly_divmod(n, d)
        for x in (0.3, -1.2, 2.5):
            lhs = poly_eval(n, x)
            rhs = poly_eval(q, x) * poly_eval(d, x) + poly_eval(r, x)
            assert abs(lhs - rhs) < 1e-6 * max(1.0, abs(lhs))


def test_division_by_constant_terminates():
    # Regression: float round-off in the leading-term cancellation once looped forever.
    q, r = poly_divmod([0.4, 0.5, -1.03, 3.6, -2.68], [1.7141])
    assert r == [0.0]
    assert len(q) == 5


def test_derivative_and_integral_inverse():
    c = [3, 1, 4, 1, 5]
    assert poly_derivative(poly_integral(c)) == [float(x) for x in c]
    assert poly_derivative([0, 0, 0, 1]) == [0, 0, 3]


def test_eval_horner():
    assert poly_eval([1, 2, 3], 2) == 17.0     # 1 + 4 + 12


def test_gcd_detects_repeated_root():
    p = [-2, 5, -4, 1]                          # (x-1)^2 (x-2)
    g = poly_gcd(p, poly_derivative(p))
    assert abs(g[0] + 1) < 1e-6 and abs(g[1] - 1) < 1e-6   # monic x - 1


def test_gcd_coprime_is_constant():
    assert poly_gcd([1, 1], [1, 0, 1]) == [1.0]


def test_validation():
    with pytest.raises(ValueError):
        poly_divmod([1, 2, 3], [0])
