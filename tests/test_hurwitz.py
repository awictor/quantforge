import math

import pytest

from quantforge import hurwitz_zeta, polygamma, riemann_zeta, digamma


def close(a, b, tol=1e-8):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_hurwitz_reduces_to_riemann():
    for s in (2.0, 3.0, 4.5, 7.0):
        assert close(hurwitz_zeta(s, 1.0), riemann_zeta(s))


def test_hurwitz_half_and_shift_identities():
    # zeta(s, 1/2) = (2^s - 1) zeta(s)
    for s in (2.0, 3.0, 5.0):
        assert close(hurwitz_zeta(s, 0.5), (2 ** s - 1) * riemann_zeta(s))
    # zeta(s, 2) = zeta(s) - 1
    for s in (2.0, 3.0, 4.0):
        assert close(hurwitz_zeta(s, 2.0), riemann_zeta(s) - 1.0)


def test_trigamma_closed_forms():
    assert close(polygamma(1, 1.0), math.pi ** 2 / 6)        # trigamma(1) = zeta(2)
    assert close(polygamma(1, 2.0), math.pi ** 2 / 6 - 1.0)
    assert close(polygamma(1, 0.5), math.pi ** 2 / 2)


def test_polygamma_higher_order():
    assert close(polygamma(2, 1.0), -2 * riemann_zeta(3))    # psi''(1) = -2 zeta(3)
    assert close(polygamma(0, 3.3), digamma(3.3))            # order 0 = digamma


def test_polygamma_matches_finite_difference():
    h = 1e-5
    for x in (1.7, 4.2):
        fd = (digamma(x + h) - digamma(x - h)) / (2 * h)
        assert close(polygamma(1, x), fd, 1e-5)
    x = 2.5
    fd2 = (digamma(x + h) - 2 * digamma(x) + digamma(x - h)) / (h * h)
    assert close(polygamma(2, x), fd2, 1e-3)


def test_domain_errors():
    with pytest.raises(ValueError):
        hurwitz_zeta(1.0, 1.0)
    with pytest.raises(ValueError):
        hurwitz_zeta(0.5, 1.0)
    with pytest.raises(ValueError):
        hurwitz_zeta(2.0, -1.0)
    with pytest.raises(ValueError):
        polygamma(-1, 1.0)
    with pytest.raises(ValueError):
        polygamma(1, 0.0)
