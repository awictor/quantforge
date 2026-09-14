import math

import pytest

from quantforge import hermite_h, hermite_he, laguerre_l, chebyshev_t, chebyshev_u


def close(a, b, tol=1e-8):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _integ(f, a, b, N=200000):
    h = (b - a) / N
    s = 0.5 * (f(a) + f(b))
    for i in range(1, N):
        s += f(a + i * h)
    return s * h


def test_hermite_physicists_explicit():
    for x in (-1.3, -0.4, 0.6, 1.7):
        assert close(hermite_h(0, x), 1.0)
        assert close(hermite_h(1, x), 2 * x)
        assert close(hermite_h(2, x), 4 * x * x - 2)
        assert close(hermite_h(3, x), 8 * x ** 3 - 12 * x)
        assert close(hermite_h(4, x), 16 * x ** 4 - 48 * x * x + 12)


def test_hermite_probabilists_and_relation():
    for x in (-1.3, 0.6, 1.7):
        assert close(hermite_he(2, x), x * x - 1)
        assert close(hermite_he(3, x), x ** 3 - 3 * x)
        for n in range(6):
            assert close(hermite_he(n, x), 2 ** (-n / 2) * hermite_h(n, x / math.sqrt(2)), 1e-7)


def test_laguerre_explicit():
    for x in (0.2, 1.5, 3.0):
        assert close(laguerre_l(0, x), 1.0)
        assert close(laguerre_l(1, x), 1 - x)
        assert close(laguerre_l(2, x), (x * x - 4 * x + 2) / 2)
        assert close(laguerre_l(3, x), (-x ** 3 + 9 * x * x - 18 * x + 6) / 6)


def test_generalized_laguerre():
    for x in (0.3, 2.0):
        assert close(laguerre_l(1, x, 1.0), 2 - x)
        assert close(laguerre_l(2, x, 1.0), (x * x - 6 * x + 6) / 2)


def test_chebyshev_trig_identities():
    for th in (0.3, 1.1, 2.4, 2.9):
        x = math.cos(th)
        for n in range(8):
            assert close(chebyshev_t(n, x), math.cos(n * th))
            assert close(chebyshev_u(n, x), math.sin((n + 1) * th) / math.sin(th), 1e-7)


def test_chebyshev_explicit():
    for x in (-0.7, 0.5):
        assert close(chebyshev_t(2, x), 2 * x * x - 1)
        assert close(chebyshev_t(3, x), 4 * x ** 3 - 3 * x)
        assert close(chebyshev_u(2, x), 4 * x * x - 1)


def test_chebyshev_orthogonality_theta():
    # int_-1^1 T_m T_n / sqrt(1-x^2) dx = int_0^pi cos(m th) cos(n th) dth
    for m in range(4):
        for n in range(4):
            val = _integ(lambda th: chebyshev_t(m, math.cos(th)) * chebyshev_t(n, math.cos(th)),
                         0.0, math.pi)
            if m != n:
                expect = 0.0
            elif m == 0:
                expect = math.pi
            else:
                expect = math.pi / 2
            assert close(val, expect, 1e-4)


@pytest.mark.slow
def test_hermite_orthogonality():
    for m in range(4):
        for n in range(4):
            val = _integ(lambda x: hermite_h(m, x) * hermite_h(n, x) * math.exp(-x * x),
                         -8.0, 8.0, N=400000)
            expect = math.sqrt(math.pi) * 2 ** n * math.factorial(n) if m == n else 0.0
            assert close(val, expect, 1e-5)


@pytest.mark.slow
def test_laguerre_orthogonality():
    for m in range(4):
        for n in range(4):
            val = _integ(lambda x: laguerre_l(m, x) * laguerre_l(n, x) * math.exp(-x),
                         0.0, 40.0, N=400000)
            expect = 1.0 if m == n else 0.0
            assert close(val, expect, 1e-4)


def test_domain_errors():
    for fn in (hermite_h, hermite_he, laguerre_l, chebyshev_t, chebyshev_u):
        with pytest.raises(ValueError):
            fn(-1, 0.5)
