"""Dated bond price / yield from calendar cashflows (bondmath module)."""

import pytest

from quantforge import (
    dated_bond_cashflows, dated_bond_price, dated_bond_yield,
    bond_cashflows, bond_price_from_yield,
)


CF = dated_bond_cashflows((2024, 1, 15), 2.0, 100, 0.05, 2, "30/360")


def test_matches_uniform_30_360():
    p = dated_bond_price((2024, 1, 15), CF, 0.04, "30/360")
    u = bond_price_from_yield(bond_cashflows(100, 0.05, 2, 2), 0.04)
    assert p == pytest.approx(u, abs=1e-6)


def test_yield_round_trip():
    p = dated_bond_price((2024, 1, 15), CF, 0.04, "30/360")
    assert dated_bond_yield((2024, 1, 15), CF, p, "30/360") == pytest.approx(
        0.04, abs=1e-8)


def test_settle_after_coupon_drops_it():
    full = dated_bond_price((2024, 1, 15), CF, 0.04, "30/360")
    later = dated_bond_price((2024, 8, 1), CF, 0.04, "30/360")
    assert later < full


def test_higher_yield_lower_price():
    assert dated_bond_price((2024, 1, 15), CF, 0.06, "30/360") < \
        dated_bond_price((2024, 1, 15), CF, 0.02, "30/360")


def test_validation():
    with pytest.raises(ValueError):
        dated_bond_yield((2024, 1, 15), CF, -1.0, "30/360")
