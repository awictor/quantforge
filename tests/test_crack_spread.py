"""Refinery crack-spread option."""

import math

import pytest

from quantforge import crack_spread_option as crack
from quantforge.commodity import bachelier_spread_option as bsp


def test_single_product_reduces_to_bachelier_spread():
    c = crack(80, [95], [1.0], 10, 5.0, [6.0], [0.6], 0.03, 1.0)
    b = bsp(95, 80, 10, 6.0, 5.0, 0.6, 0.03, 1.0)
    assert abs(c - b) < 1e-9


def test_321_crack_positive():
    cr = crack(80, [95, 90], [2/3, 1/3], 5, 5.0, [7.0, 6.5], [0.7, 0.65], 0.03, 1.0)
    assert cr > 0


def test_put_call_parity():
    args = (80, [95, 90], [2/3, 1/3], 5, 5.0, [7.0, 6.5], [0.7, 0.65], 0.03, 1.0)
    c = crack(*args, is_call=True)
    p = crack(*args, is_call=False)
    basket = 2/3 * 95 + 1/3 * 90
    assert abs((c - p) - math.exp(-0.03) * (basket - 80 - 5)) < 1e-9


def test_higher_product_vol_raises():
    base = crack(80, [95, 90], [2/3, 1/3], 5, 5.0, [7.0, 6.5], [0.7, 0.65], 0.03, 1.0)
    high = crack(80, [95, 90], [2/3, 1/3], 5, 5.0, [10.0, 9.0], [0.7, 0.65], 0.03, 1.0)
    assert high > base


def test_validation():
    with pytest.raises(ValueError):
        crack(80, [95], [1, 2], 5, 5.0, [6.0], [0.6], 0.03, 1.0)   # length mismatch
