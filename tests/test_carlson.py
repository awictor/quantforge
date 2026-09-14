import math

import pytest

from quantforge import (
    carlson_rf,
    carlson_rc,
    carlson_rd,
    carlson_rj,
    elliptic_f,
    elliptic_e_incomplete,
    elliptic_pi,
    elliptic_k,
    elliptic_e,
)


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_symmetric_form_spot_values():
    assert close(carlson_rf(2, 2, 2), 1 / math.sqrt(2))
    assert close(carlson_rf(0, 1, 1), math.pi / 2)      # = RC(0, 1)
    assert close(carlson_rc(3, 3), 1 / math.sqrt(3))
    assert close(carlson_rc(0, 4), math.pi / 4)         # pi/(2 sqrt y)
    x, y = 0.3, 1.7
    assert close(carlson_rc(x, y), math.acos(math.sqrt(x / y)) / math.sqrt(y - x))
    assert close(carlson_rd(2, 2, 2), 2 ** -1.5)
    assert close(carlson_rj(2, 2, 2, 2), 2 ** -1.5)


def test_rj_reduces_to_rd():
    # R_J(x, y, z, z) = R_D(x, y, z)
    assert close(carlson_rj(0.5, 1.0, 2.0, 2.0), carlson_rd(0.5, 1.0, 2.0))


def test_rf_symmetry():
    assert close(carlson_rf(1, 2, 3), carlson_rf(3, 1, 2))
    assert close(carlson_rf(1, 2, 3), carlson_rf(2, 3, 1))


def _F_int(phi, m, N=200000):
    h = phi / N

    def f(t):
        return 1.0 / math.sqrt(1.0 - m * math.sin(t) ** 2)

    s = 0.5 * (f(0.0) + f(phi))
    for i in range(1, N):
        s += f(i * h)
    return s * h


def _E_int(phi, m, N=200000):
    h = phi / N

    def f(t):
        return math.sqrt(1.0 - m * math.sin(t) ** 2)

    s = 0.5 * (f(0.0) + f(phi))
    for i in range(1, N):
        s += f(i * h)
    return s * h


def test_incomplete_f_and_e_vs_integration():
    for phi, m in ((0.7, 0.3), (1.2, 0.6), (1.0, 0.9), (1.4, 0.5)):
        assert close(elliptic_f(phi, m), _F_int(phi, m), 1e-6)
        assert close(elliptic_e_incomplete(phi, m), _E_int(phi, m), 1e-6)


def test_incomplete_reduces_to_complete():
    for m in (0.1, 0.5, 0.8):
        assert close(elliptic_f(math.pi / 2, m), elliptic_k(m))
        assert close(elliptic_e_incomplete(math.pi / 2, m), elliptic_e(m))


def _Pi_int(n, phi, m, N=200000):
    h = phi / N

    def f(t):
        s2 = math.sin(t) ** 2
        return 1.0 / ((1.0 - n * s2) * math.sqrt(1.0 - m * s2))

    s = 0.5 * (f(0.0) + f(phi))
    for i in range(1, N):
        s += f(i * h)
    return s * h


def test_third_kind_vs_integration():
    for n, phi, m in ((0.3, 0.8, 0.4), (0.5, 1.1, 0.6), (-0.4, 1.0, 0.5), (0.2, 1.3, 0.2)):
        assert close(elliptic_pi(n, phi, m), _Pi_int(n, phi, m), 1e-6)


def test_third_kind_zero_char_is_first_kind():
    assert close(elliptic_pi(0.0, 1.0, 0.5), elliptic_f(1.0, 0.5))


def test_domain_errors():
    with pytest.raises(ValueError):
        carlson_rf(-1, 1, 1)
    with pytest.raises(ValueError):
        carlson_rf(0, 0, 1)     # two zeros
    with pytest.raises(ValueError):
        carlson_rd(1, 1, 0)     # z must be positive
    with pytest.raises(ValueError):
        carlson_rj(1, 1, 1, -1)  # p > 0 only
