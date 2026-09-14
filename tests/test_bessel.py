"""Tests for Bessel functions J0/J1/Jn/Y0/Y1, vs known values, zeros, and recurrence."""

import math
import random

import pytest

from quantforge.bessel import bessel_j0, bessel_j1, bessel_jn, bessel_y0, bessel_y1


def _close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1, abs(a), abs(b))


def _jn_integral(n, x, N=100000):
    h = math.pi / N
    s = 0.0
    for i in range(N + 1):
        t = i * h
        w = 0.5 if i in (0, N) else 1.0
        s += w * math.cos(n * t - x * math.sin(t))
    return s * h / math.pi


def test_j_known_values():
    assert _close(bessel_j0(0), 1.0)
    assert _close(bessel_j1(0), 0.0, 1e-9)
    assert _close(bessel_j0(1), 0.7651976865579666)
    assert _close(bessel_j1(1), 0.4400505857449335)
    assert _close(bessel_j0(5), -0.1775967713143383)
    assert _close(bessel_j1(5), -0.3275791375914652)
    assert _close(bessel_j0(10), -0.2459357644513483)


def test_y_known_values():
    assert _close(bessel_y0(1), 0.0882569642156769)
    assert _close(bessel_y1(1), -0.7812128213002887)
    assert _close(bessel_y0(5), -0.3085176252490338)


def test_j_zeros():
    assert abs(bessel_j0(2.404825557695773)) < 1e-6
    assert abs(bessel_j0(5.520078110286311)) < 1e-6
    assert abs(bessel_j1(3.831705970207512)) < 1e-6


def test_fuzz_recurrence():
    rng = random.Random(561)
    for _ in range(3000):
        x = rng.uniform(0.5, 30)
        n = rng.randint(1, 10)
        lhs = bessel_jn(n - 1, x) + bessel_jn(n + 1, x)
        rhs = 2 * n / x * bessel_jn(n, x)
        assert _close(lhs, rhs, 1e-7)


def test_jn_vs_integral():
    for n, x in [(2, 3), (3, 5), (4, 10), (6, 15), (3, 25)]:
        assert _close(bessel_jn(n, x), _jn_integral(n, x), 1e-6)


def test_jn_matches_j0_j1():
    for x in (0.5, 3, 7, 15):
        assert _close(bessel_jn(0, x), bessel_j0(x))
        assert _close(bessel_jn(1, x), bessel_j1(x))


def test_jn_specific():
    assert _close(bessel_jn(2, 1), 0.11490348493190049, 1e-8)
    assert _close(bessel_jn(3, 5), _jn_integral(3, 5), 1e-8)


def test_jn_at_zero():
    assert bessel_jn(3, 0) == 0.0
    assert bessel_jn(5, 0) == 0.0


def test_negative_argument_parity():
    # J_n(-x) = (-1)^n J_n(x)
    for n in range(5):
        assert _close(bessel_jn(n, -3.0), (-1) ** n * bessel_jn(n, 3.0), 1e-7)


def test_large_argument():
    assert _close(bessel_j0(20), 0.1670246643405831, 1e-5)
    assert _close(bessel_j1(20), 0.06683312417584991, 1e-5)


def test_domain_errors():
    with pytest.raises(ValueError):
        bessel_y0(0)
    with pytest.raises(ValueError):
        bessel_y1(-1)
    with pytest.raises(ValueError):
        bessel_jn(-1, 1)
