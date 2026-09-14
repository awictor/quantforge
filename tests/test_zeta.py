import math

import pytest

from quantforge import riemann_zeta, dirichlet_eta


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_zeta_even_closed_forms():
    assert close(riemann_zeta(2), math.pi ** 2 / 6)
    assert close(riemann_zeta(4), math.pi ** 4 / 90)
    assert close(riemann_zeta(6), math.pi ** 6 / 945)
    assert close(riemann_zeta(8), math.pi ** 8 / 9450)


def test_zeta_apery_and_large():
    assert close(riemann_zeta(3), 1.2020569031595942)
    assert close(riemann_zeta(10), 1.0009945751278181, 1e-10)
    assert close(riemann_zeta(50), 1.0, 1e-12)


def test_eta_closed_forms():
    assert close(dirichlet_eta(1), math.log(2))
    assert close(dirichlet_eta(2), math.pi ** 2 / 12)
    assert close(dirichlet_eta(4), 7 * math.pi ** 4 / 720)


def test_eta_zeta_relation():
    for s in (1.5, 2.0, 3.0, 5.0, 8.0):
        assert close(dirichlet_eta(s), (1.0 - 2.0 ** (1.0 - s)) * riemann_zeta(s))


def test_zeta_critical_strip_continuation():
    # values on (0, 1) via the eta relation (DLMF reference values)
    assert close(riemann_zeta(0.5), -1.4603545088095868, 1e-8)
    assert close(riemann_zeta(0.25), -0.8132784052618923, 1e-8)


def test_zeta_pole_and_domain():
    with pytest.raises(ValueError):
        riemann_zeta(1.0)
    with pytest.raises(ValueError):
        dirichlet_eta(0.0)
    with pytest.raises(ValueError):
        dirichlet_eta(-1.0)


@pytest.mark.slow
def test_zeta_matches_direct_sum():
    for s in (2.5, 3.7):
        direct = sum(n ** (-s) for n in range(1, 2_000_000))
        assert close(riemann_zeta(s), direct, 1e-5)
