"""Black-Karasinski short-rate model."""

import math

import pytest

from quantforge import bk_zero_coupon_bond as bk


def test_bond_in_unit_interval():
    p = bk(0.03, 0.1, math.log(0.03), 0.2, 5.0, 80)
    assert 0.0 < p < 1.0


def test_order_of_magnitude_near_flat_discount():
    p = bk(0.03, 0.1, math.log(0.03), 0.2, 5.0, 80)
    assert abs(p - math.exp(-0.03 * 5)) < 0.1


def test_higher_rate_lower_price():
    lo = bk(0.03, 0.1, math.log(0.03), 0.2, 5.0, 80)
    hi = bk(0.06, 0.1, math.log(0.03), 0.2, 5.0, 80)
    assert hi < lo


def test_short_maturity_near_one():
    assert bk(0.03, 0.1, math.log(0.03), 0.2, 0.01, 20) > 0.999


def test_converges_across_steps():
    a = bk(0.03, 0.1, math.log(0.03), 0.2, 5.0, 40)
    b = bk(0.03, 0.1, math.log(0.03), 0.2, 5.0, 160)
    assert abs(a - b) < 0.02


def test_higher_vol_lowers_price():
    # Log-normal: higher sigma raises E[r] = E[exp(x)] (Jensen), lowering the bond.
    lo_vol = bk(0.03, 0.1, math.log(0.03), 0.2, 5.0, 80)
    hi_vol = bk(0.03, 0.1, math.log(0.03), 0.5, 5.0, 80)
    assert hi_vol < lo_vol


def test_validation():
    with pytest.raises(ValueError):
        bk(-0.01, 0.1, 0.0, 0.2, 5.0)
    with pytest.raises(ValueError):
        bk(0.03, 0.1, 0.0, 0.2, 0.0)
