import math
import random

import pytest

from quantforge import jacobi_sn, jacobi_cn, jacobi_dn, jacobi_am, elliptic_k


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_pythagorean_identities():
    random.seed(1)
    for _ in range(500):
        m = random.random() * 0.999
        u = random.uniform(-8, 8)
        sn, cn, dn = jacobi_sn(u, m), jacobi_cn(u, m), jacobi_dn(u, m)
        assert close(sn * sn + cn * cn, 1.0)
        assert close(dn * dn + m * sn * sn, 1.0)


def test_m_zero_is_trig():
    for u in (0.3, 1.7, -2.2, 5.0):
        assert close(jacobi_sn(u, 0.0), math.sin(u))
        assert close(jacobi_cn(u, 0.0), math.cos(u))
        assert close(jacobi_dn(u, 0.0), 1.0)


def test_m_one_is_hyperbolic():
    # sn(u,1)=tanh u, cn(u,1)=dn(u,1)=sech u
    for u in (0.5, 1.5, -1.0):
        assert close(jacobi_sn(u, 1.0), math.tanh(u))
        assert close(jacobi_cn(u, 1.0), 1.0 / math.cosh(u))
        assert close(jacobi_dn(u, 1.0), 1.0 / math.cosh(u))


def test_quarter_period_values():
    for m in (0.1, 0.5, 0.8, 0.99):
        K = elliptic_k(m)
        assert close(jacobi_sn(K, m), 1.0, 1e-8)
        assert close(jacobi_cn(K, m), 0.0, 1e-7)
        assert close(jacobi_dn(K, m), math.sqrt(1.0 - m), 1e-8)


def test_amplitude_inverts_incomplete_integral():
    def F(phi, m, N=100000):
        h = phi / N

        def f(t):
            return 1.0 / math.sqrt(1.0 - m * math.sin(t) ** 2)

        s = 0.5 * (f(0.0) + f(phi))
        for i in range(1, N):
            s += f(i * h)
        return s * h

    for u, m in ((0.7, 0.3), (1.2, 0.6), (0.4, 0.9), (2.0, 0.5)):
        K = elliptic_k(m)
        if abs(u) < K:
            assert close(F(jacobi_am(u, m), m), u, 1e-4)


def test_derivative_sn_is_cn_dn():
    h = 1e-6
    for u, m in ((0.8, 0.4), (1.5, 0.7)):
        d = (jacobi_sn(u + h, m) - jacobi_sn(u - h, m)) / (2 * h)
        assert close(d, jacobi_cn(u, m) * jacobi_dn(u, m), 1e-5)


def test_period_four_k():
    for m in (0.3, 0.7):
        K = elliptic_k(m)
        for u in (0.5, 1.3):
            assert close(jacobi_sn(u + 4 * K, m), jacobi_sn(u, m), 1e-7)


def test_domain_errors():
    with pytest.raises(ValueError):
        jacobi_am(0.5, -0.1)
    with pytest.raises(ValueError):
        jacobi_am(0.5, 1.1)
