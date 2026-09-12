"""Bivariate copulas: Gaussian, Clayton, Gumbel."""

import pytest

from quantforge import (
    gaussian_copula, clayton_copula, gumbel_copula,
    clayton_lower_tail_dependence, gumbel_upper_tail_dependence,
    clayton_theta_from_tau, gumbel_theta_from_tau,
)


def test_boundary_zero():
    assert gaussian_copula(0.5, 0, 0.5) == 0
    assert clayton_copula(0.5, 0, 2) == 0
    assert gumbel_copula(0.5, 0, 2) == 0


def test_boundary_uniform_margin():
    assert gaussian_copula(0.5, 1, 0.5) == pytest.approx(0.5)
    assert clayton_copula(0.5, 1, 2) == pytest.approx(0.5)
    assert gumbel_copula(0.5, 1, 2) == pytest.approx(0.5)


def test_independence_limits():
    assert gaussian_copula(0.5, 0.6, 0.0) == pytest.approx(0.3, abs=1e-9)
    assert clayton_copula(0.5, 0.6, 1e-8) == pytest.approx(0.3, abs=1e-6)
    assert gumbel_copula(0.5, 0.6, 1.0) == pytest.approx(0.3, abs=1e-12)


def test_frechet_bounds_and_positive_dependence():
    u, v = 0.4, 0.7
    lo, hi = max(u + v - 1, 0), min(u, v)
    for val in (gaussian_copula(u, v, 0.5), clayton_copula(u, v, 2),
                gumbel_copula(u, v, 2)):
        assert lo - 1e-9 <= val <= hi + 1e-9
        assert val > u * v   # positive dependence


def test_clayton_lower_tail_dependence():
    assert clayton_lower_tail_dependence(2) == pytest.approx(2 ** (-0.5))
    assert 0 < clayton_lower_tail_dependence(2) < 1


def test_gumbel_upper_tail_dependence():
    assert gumbel_upper_tail_dependence(2) == pytest.approx(2 - 2 ** 0.5)
    assert gumbel_upper_tail_dependence(1.0) == pytest.approx(0.0)


def test_theta_from_tau_round_trips():
    th = clayton_theta_from_tau(0.5)
    assert th / (th + 2) == pytest.approx(0.5)
    thg = gumbel_theta_from_tau(0.5)
    assert 1 - 1 / thg == pytest.approx(0.5)


def test_validation():
    with pytest.raises(ValueError):
        clayton_copula(0.5, 0.5, -1)
    with pytest.raises(ValueError):
        gumbel_copula(0.5, 0.5, 0.5)
    with pytest.raises(ValueError):
        clayton_theta_from_tau(1.0)
