"""Dated accrued interest and clean price (bondmath module)."""

import pytest

from quantforge import (
    dated_bond_cashflows, dated_bond_price, dated_clean_price,
    dated_accrued_interest,
)


CF = dated_bond_cashflows((2024, 1, 15), 2.0, 100, 0.05, 2, "30/360")
PREV, NEXT = (2024, 1, 15), (2024, 7, 15)


def test_accrued_half_period():
    # 3 months into a 6-month period: half of the 2.5 coupon.
    assert dated_accrued_interest((2024, 4, 15), PREV, NEXT, 100, 0.05, 2,
                                  "30/360") == pytest.approx(1.25, abs=1e-9)


def test_accrued_zero_at_coupon():
    assert dated_accrued_interest(PREV, PREV, NEXT, 100, 0.05, 2, "30/360") == 0.0


def test_clean_is_dirty_minus_accrued():
    settle = (2024, 4, 15)
    dirty = dated_bond_price(settle, CF, 0.04, "30/360")
    accrued = dated_accrued_interest(settle, PREV, NEXT, 100, 0.05, 2, "30/360")
    clean = dated_clean_price(settle, CF, 0.04, PREV, NEXT, 100, 0.05, 2, "30/360")
    assert dirty - clean == pytest.approx(accrued, abs=1e-9)


def test_clean_equals_dirty_at_coupon():
    dirty = dated_bond_price(PREV, CF, 0.04, "30/360")
    clean = dated_clean_price(PREV, CF, 0.04, PREV, NEXT, 100, 0.05, 2, "30/360")
    assert clean == pytest.approx(dirty, abs=1e-9)


def test_accrued_grows_toward_next_coupon():
    early = dated_accrued_interest((2024, 2, 15), PREV, NEXT, 100, 0.05, 2, "30/360")
    late = dated_accrued_interest((2024, 6, 15), PREV, NEXT, 100, 0.05, 2, "30/360")
    assert late > early


def test_validation():
    with pytest.raises(ValueError):
        dated_accrued_interest((2023, 12, 1), PREV, NEXT, 100, 0.05, 2)
