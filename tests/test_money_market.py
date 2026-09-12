"""Money-market yield conventions."""

import pytest

from quantforge import (
    price_from_discount, bank_discount_yield, money_market_yield,
    bond_equivalent_yield, discount_to_bond_equivalent, holding_period_return,
)


FACE, D, DAYS = 100.0, 0.05, 90


def test_price_discount_round_trip():
    p = price_from_discount(FACE, D, DAYS)
    assert p == pytest.approx(98.75)
    assert bank_discount_yield(FACE, p, DAYS) == pytest.approx(D)


def test_yield_ordering():
    p = price_from_discount(FACE, D, DAYS)
    disc = bank_discount_yield(FACE, p, DAYS)
    mm = money_market_yield(FACE, p, DAYS)
    bey = bond_equivalent_yield(FACE, p, DAYS)
    assert mm > disc          # divides by price, not face
    assert bey > mm           # 365 vs 360
    assert bey > D


def test_discount_to_bond_equivalent_above_input():
    assert discount_to_bond_equivalent(D, DAYS) > D


def test_holding_period_return():
    assert holding_period_return(98, 100, 1) == pytest.approx(3 / 98)


def test_longer_maturity_lower_price():
    assert price_from_discount(FACE, D, 180) < price_from_discount(FACE, D, 90)


def test_validation():
    with pytest.raises(ValueError):
        price_from_discount(FACE, D, 0)
    with pytest.raises(ValueError):
        holding_period_return(0, 100)
