import math

import pytest

from quantforge import polylog, dilog, riemann_zeta, dirichlet_eta


def close(a, b, tol=1e-10):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _brute_dilog(x, n=200000):
    s = 0.0
    xk = x
    for k in range(1, n):
        s += xk / (k * k)
        xk *= x
    return s


def test_dilog_closed_forms():
    assert close(dilog(1.0), math.pi ** 2 / 6)
    assert close(dilog(-1.0), -math.pi ** 2 / 12)
    assert close(dilog(0.5), math.pi ** 2 / 12 - math.log(2) ** 2 / 2)
    assert close(dilog(0.0), 0.0)


def test_dilog_matches_series_in_disk():
    for x in (0.9, -0.9, 0.3, -0.5, 0.7):
        assert close(dilog(x), _brute_dilog(x), 1e-8)


def test_dilog_inversion_branch():
    # for x < -1 the inversion identity must match numerical integration of -ln(1-t)/t
    def integ(x, N=1_000_000):
        h = x / N

        def f(t):
            return 1.0 if t == 0.0 else -math.log(1.0 - t) / t

        s = 0.5 * (f(0.0) + f(x))
        for i in range(1, N):
            s += f(i * h)
        return s * h

    for x in (-2.0, -5.0):
        assert close(dilog(x), integ(x), 1e-4)


def test_polylog_endpoints():
    assert close(polylog(2, 1.0), riemann_zeta(2))
    assert close(polylog(3, 1.0), riemann_zeta(3))
    assert close(polylog(2, -1.0), -dirichlet_eta(2))
    assert close(polylog(2, 0.5), dilog(0.5))


def test_polylog_order_one_is_log():
    for z in (0.3, -0.5, 0.8, -0.9):
        assert close(polylog(1, z), -math.log(1.0 - z), 1e-9)


def test_domain_errors():
    with pytest.raises(ValueError):
        polylog(2, 1.5)
    with pytest.raises(ValueError):
        polylog(2, -1.5)
    with pytest.raises(ValueError):
        dilog(2.0)
    with pytest.raises(ValueError):
        polylog(1, 1.0)     # Li_1(1) diverges -> zeta(1) pole
