"""Schedule-driven bond cashflows with day-count accruals (bondmath module)."""

import pytest

from quantforge import dated_bond_cashflows


def test_semiannual_30_360_structure():
    cf = dated_bond_cashflows((2024, 1, 15), 2.0, 100, 0.05, 2, "30/360")
    assert len(cf) == 4
    # every 30/360 semiannual period is exactly half a year
    for _, tau, _ in cf:
        assert tau == pytest.approx(0.5, abs=1e-9)
    # coupon = face * rate * tau = 2.5
    assert cf[0][2] == pytest.approx(2.5, abs=1e-9)


def test_final_flow_includes_face():
    cf = dated_bond_cashflows((2024, 1, 15), 2.0, 100, 0.05, 2, "30/360")
    assert cf[-1][2] > 100.0
    assert cf[-1][0] == (2026, 1, 15)


def test_act360_accrual_exceeds_half():
    # act/360 over ~182 days is slightly more than 0.5.
    cf = dated_bond_cashflows((2024, 1, 15), 1.0, 100, 0.05, 2, "act/360")
    assert cf[0][1] > 0.5


def test_zero_coupon_only_face():
    cf = dated_bond_cashflows((2024, 1, 15), 1.0, 100, 0.0, 1, "30/360")
    assert cf[-1][2] == pytest.approx(100.0, abs=1e-9)


def test_end_of_month_roll():
    cf = dated_bond_cashflows((2024, 1, 31), 1.0, 100, 0.05, 2, "30/360",
                              end_of_month=True)
    # month-end schedule dates
    assert cf[0][0] == (2024, 7, 31)


def test_validation():
    with pytest.raises(ValueError):
        dated_bond_cashflows((2024, 1, 15), 1.0, 100, -0.01)
