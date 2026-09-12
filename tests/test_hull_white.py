"""Hull-White (extended Vasicek) short-rate model."""

import math

import pytest

from quantforge import hw_zero_from_curve as hw, hw_B


def test_refits_flat_curve_at_time_zero():
    r = 0.03
    P0 = lambda T: math.exp(-r * T)
    for T in (1, 5, 10, 20):
        assert abs(hw(P0, r, 0.1, 0.01, 0.0, T) - P0(T)) < 1e-4


def test_refits_upward_curve():
    P0 = lambda T: math.exp(-(0.02 + 0.005 * T) * T) if T > 0 else 1.0
    # initial instantaneous forward at 0 ~ 0.02
    assert abs(hw(P0, 0.02, 0.1, 0.01, 0.0, 7.0) - P0(7.0)) < 0.01


def test_ho_lee_limit_refits():
    r = 0.03
    P0 = lambda T: math.exp(-r * T)
    assert abs(hw(P0, r, 1e-9, 0.01, 0.0, 5.0) - P0(5.0)) < 1e-3


def test_B_factor_limit_and_sign():
    assert abs(hw_B(1e-9, 5.0) - 5.0) < 1e-3     # -> tau as a -> 0
    assert hw_B(0.1, 5.0) > 0


def test_discount_in_unit_interval():
    r = 0.03
    P0 = lambda T: math.exp(-r * T)
    assert 0.0 < hw(P0, r, 0.1, 0.01, 2.0, 10.0) < 1.0


def test_validation():
    P0 = lambda T: math.exp(-0.03 * T)
    with pytest.raises(ValueError):
        hw(P0, 0.03, 0.1, -0.01, 0.0, 5.0)       # negative sigma
    with pytest.raises(ValueError):
        hw(P0, 0.03, 0.1, 0.01, 5.0, 2.0)        # T < t
