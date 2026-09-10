"""Tests for Lee moment-bound wing checks on an SVI slice."""

import pytest

from quantforge import SVIParams, lee_wing_slopes, lee_bounds_ok


def test_slopes_are_b_times_one_plus_minus_rho():
    p = SVIParams(a=0.04, b=0.5, rho=-0.3, m=0.0, s=0.1)
    left, right = lee_wing_slopes(p)
    assert left == pytest.approx(0.5 * (1 - (-0.3)))
    assert right == pytest.approx(0.5 * (1 + (-0.3)))


def test_slopes_match_asymptotic_total_variance():
    # At large |k| the SVI total variance is linear with slope b(1 +/- rho).
    p = SVIParams(a=0.04, b=0.5, rho=-0.3, m=0.0, s=0.1)
    _, right = lee_wing_slopes(p)
    emp = p.total_variance(50.0) - p.total_variance(49.0)
    assert emp == pytest.approx(right, abs=1e-3)


def test_reasonable_slice_passes():
    p = SVIParams(a=0.04, b=0.5, rho=-0.3, m=0.0, s=0.1)
    assert lee_bounds_ok(p)


def test_steep_wing_fails():
    p = SVIParams(a=0.04, b=2.0, rho=-0.5, m=0.0, s=0.1)  # right slope = 3 > 2
    assert not lee_bounds_ok(p)


def test_consistent_with_is_arbitrage_free_wings():
    for b, rho in ((0.3, -0.4), (0.8, 0.2), (1.5, -0.5), (2.0, 0.5)):
        p = SVIParams(a=0.04, b=b, rho=rho, m=0.0, s=0.1)
        assert lee_bounds_ok(p) == p.is_arbitrage_free_wings()
