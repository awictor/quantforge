"""Accrued interest and clean/dirty bond prices (bondmath module)."""

import pytest

from quantforge import (
    bond_cashflows, accrued_interest, clean_price, dirty_price,
    bond_price_from_yield,
)


CF = bond_cashflows(100.0, 0.05, 5.0, freq=2)
Y = 0.04


def test_accrued_half_period():
    # Half of a 2.5 semiannual coupon.
    assert accrued_interest(100.0, 0.05, 2, 0.5) == pytest.approx(1.25, abs=1e-12)


def test_accrued_zero_at_coupon_date():
    assert accrued_interest(100.0, 0.05, 2, 0.0) == 0.0


def test_dirty_equals_present_value():
    assert dirty_price(CF, Y) == pytest.approx(bond_price_from_yield(CF, Y), abs=1e-12)


def test_clean_is_dirty_minus_accrued():
    d = dirty_price(CF, Y)
    c = clean_price(CF, Y, 100.0, 0.05, 2, 0.5)
    assert d - c == pytest.approx(1.25, abs=1e-9)


def test_clean_equals_dirty_at_coupon():
    assert clean_price(CF, Y, 100.0, 0.05, 2, 0.0) == pytest.approx(
        dirty_price(CF, Y), abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        accrued_interest(100.0, 0.05, 2, 1.5)   # fraction out of range
    with pytest.raises(ValueError):
        accrued_interest(100.0, 0.05, 0, 0.5)   # bad freq
