import math

import pytest

from quantforge import legendre_p, assoc_legendre, spherical_harmonic_real


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _integ(f, a, b, N=100000):
    h = (b - a) / N
    s = 0.5 * (f(a) + f(b))
    for i in range(1, N):
        s += f(a + i * h)
    return s * h


def test_legendre_explicit_polynomials():
    for x in (-0.8, -0.3, 0.0, 0.4, 0.9):
        assert close(legendre_p(0, x), 1.0)
        assert close(legendre_p(1, x), x)
        assert close(legendre_p(2, x), 0.5 * (3 * x * x - 1))
        assert close(legendre_p(3, x), 0.5 * (5 * x ** 3 - 3 * x))
        assert close(legendre_p(4, x), (35 * x ** 4 - 30 * x * x + 3) / 8)
        assert close(legendre_p(5, x), (63 * x ** 5 - 70 * x ** 3 + 15 * x) / 8)


def test_legendre_endpoints():
    for l in range(8):
        assert close(legendre_p(l, 1.0), 1.0)
        assert close(legendre_p(l, -1.0), (-1) ** l)


def test_legendre_orthogonality():
    for l in range(5):
        for k in range(5):
            val = _integ(lambda x: legendre_p(l, x) * legendre_p(k, x), -1.0, 1.0)
            expect = 2.0 / (2 * l + 1) if l == k else 0.0
            assert close(val, expect, 1e-5)


def test_assoc_reduces_to_legendre():
    for l in range(6):
        for x in (-0.7, 0.2, 0.85):
            assert close(assoc_legendre(l, 0, x), legendre_p(l, x))


def test_assoc_explicit():
    for x in (-0.6, 0.0, 0.5):
        assert close(assoc_legendre(1, 1, x), -math.sqrt(1 - x * x))
    for x in (-0.6, 0.3, 0.8):
        assert close(assoc_legendre(2, 1, x), -3 * x * math.sqrt(1 - x * x))
        assert close(assoc_legendre(2, 2, x), 3 * (1 - x * x))
    for x in (-0.5, 0.4):
        assert close(assoc_legendre(3, 2, x), 15 * x * (1 - x * x))


def _sphere_dot(l1, m1, l2, m2, Nt=200, Np=200):
    s = 0.0
    for i in range(Nt):
        th = (i + 0.5) * math.pi / Nt
        w = math.sin(th) * (math.pi / Nt) * (2 * math.pi / Np)
        for j in range(Np):
            ph = (j + 0.5) * 2 * math.pi / Np
            s += (spherical_harmonic_real(l1, m1, th, ph)
                  * spherical_harmonic_real(l2, m2, th, ph) * w)
    return s


def test_y00_is_constant():
    assert close(spherical_harmonic_real(0, 0, 1.0, 2.0), 1 / math.sqrt(4 * math.pi))


@pytest.mark.slow
def test_spherical_harmonic_orthonormal():
    modes = [(0, 0), (1, -1), (1, 0), (1, 1), (2, 0), (2, 1), (2, -2)]
    for a in modes:
        for b in modes:
            v = _sphere_dot(*a, *b)
            expect = 1.0 if a == b else 0.0
            assert close(v, expect, 6e-3)


def test_orthonormal_small_grid():
    # fast coverage: a couple of pairs on a coarse grid
    assert close(_sphere_dot(1, 0, 1, 0, 120, 120), 1.0, 1e-2)
    assert abs(_sphere_dot(1, 0, 2, 0, 120, 120)) < 1e-2


def test_domain_errors():
    with pytest.raises(ValueError):
        assoc_legendre(1, 2, 0.5)
    with pytest.raises(ValueError):
        assoc_legendre(2, -3, 0.5)
    with pytest.raises(ValueError):
        assoc_legendre(1, 1, 1.5)
    with pytest.raises(ValueError):
        spherical_harmonic_real(1, 2, 0.5, 0.5)
