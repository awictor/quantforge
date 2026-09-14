import math

import pytest

from quantforge import (
    spherical_bessel_j,
    spherical_bessel_y,
    bessel_in,
    bessel_kn,
)


def close(a, b, tol=1e-6):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_spherical_j_explicit():
    for x in (0.5, 1.3, 3.0, 7.0):
        assert close(spherical_bessel_j(0, x), math.sin(x) / x)
        assert close(spherical_bessel_j(1, x), math.sin(x) / x ** 2 - math.cos(x) / x)
        j2 = (3 / x ** 2 - 1) * math.sin(x) / x - 3 * math.cos(x) / x ** 2
        assert close(spherical_bessel_j(2, x), j2)
        j3 = (15 / x ** 3 - 6 / x) * math.sin(x) / x - (15 / x ** 2 - 1) * math.cos(x) / x
        assert close(spherical_bessel_j(3, x), j3)


def test_spherical_y_explicit():
    for x in (0.5, 1.3, 3.0, 7.0):
        assert close(spherical_bessel_y(0, x), -math.cos(x) / x)
        y2 = (-3 / x ** 2 + 1) * math.cos(x) / x - 3 * math.sin(x) / x ** 2
        assert close(spherical_bessel_y(2, x), y2)


def test_spherical_j_at_zero_and_high_order():
    assert spherical_bessel_j(0, 0.0) == 1.0
    assert spherical_bessel_j(5, 0.0) == 0.0
    v = spherical_bessel_j(10, 1.0)
    assert 0.0 < v < 1e-8


def test_spherical_wronskian():
    # j_n(x) y_n'(x) - j_n'(x) y_n(x) = 1/x^2
    def jd(n, x):
        return (n / x) * spherical_bessel_j(n, x) - spherical_bessel_j(n + 1, x)

    def yd(n, x):
        return (n / x) * spherical_bessel_y(n, x) - spherical_bessel_y(n + 1, x)

    for x in (0.8, 2.5, 5.0):
        for n in (0, 1, 3, 5):
            w = spherical_bessel_j(n, x) * yd(n, x) - jd(n, x) * spherical_bessel_y(n, x)
            assert close(w, 1 / x ** 2)


def _i_series(n, x, K=90):
    s = 0.0
    for k in range(K):
        s += (x / 2) ** (2 * k + n) / (math.factorial(k) * math.factorial(k + n))
    return s


def test_modified_i_vs_series():
    for x in (0.5, 2.0, 5.0, 12.0):
        for n in (0, 1, 2, 3, 5, 8):
            assert close(bessel_in(n, x), _i_series(n, x))


def test_modified_i_parity():
    assert close(bessel_in(2, -1.5), bessel_in(2, 1.5))
    assert close(bessel_in(3, -1.5), -bessel_in(3, 1.5))
    assert bessel_in(0, 0.0) == 1.0
    assert bessel_in(3, 0.0) == 0.0


def _k_int(n, x, N=100000, T=30.0):
    h = T / N

    def f(t):
        return math.exp(-x * math.cosh(t)) * math.cosh(n * t)

    s = 0.5 * (f(0.0) + f(T))
    for i in range(1, N):
        s += f(i * h)
    return s * h


def test_modified_k_vs_integral():
    for x in (0.5, 1.0, 2.0, 4.0):
        for n in (0, 1, 2, 3):
            assert close(bessel_kn(n, x), _k_int(n, x), 1e-4)


def test_modified_wronskian():
    # I_n(x) K_n'(x) - I_n'(x) K_n(x) = -1/x
    def Id(n, x):
        lo = bessel_in(n - 1, x) if n >= 1 else bessel_in(1, x)
        return 0.5 * (lo + bessel_in(n + 1, x))

    def Kd(n, x):
        lo = bessel_kn(n - 1, x) if n >= 1 else bessel_kn(1, x)
        return -0.5 * (lo + bessel_kn(n + 1, x))

    for x in (0.7, 1.5, 3.0):
        for n in (0, 1, 2, 4):
            w = bessel_in(n, x) * Kd(n, x) - Id(n, x) * bessel_kn(n, x)
            assert close(w, -1 / x, 1e-5)


def test_domain_errors():
    with pytest.raises(ValueError):
        bessel_in(-1, 0.5)
    with pytest.raises(ValueError):
        bessel_kn(1, 0.0)
    with pytest.raises(ValueError):
        spherical_bessel_y(1, 0.0)
    with pytest.raises(ValueError):
        spherical_bessel_j(-1, 1.0)
