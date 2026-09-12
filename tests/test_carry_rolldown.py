"""Bond carry and roll-down return decomposition."""

import pytest

from quantforge import (
    carry_return, rolldown_return, total_carry_rolldown,
    SplineZeroCurve, bond_cashflows,
)


CF = bond_cashflows(face=100, coupon_rate=0.04, maturity=5, freq=2)
UP = SplineZeroCurve([0.5, 1, 2, 3, 5, 10], [0.02, 0.025, 0.03, 0.033, 0.037, 0.04])
FLAT = SplineZeroCurve([0.5, 1, 2, 3, 5, 10], [0.03] * 6)
STEEP = SplineZeroCurve([0.5, 1, 2, 3, 5, 10], [0.01, 0.02, 0.03, 0.04, 0.05, 0.06])
INV = SplineZeroCurve([0.5, 1, 2, 3, 5, 10], [0.06, 0.05, 0.04, 0.035, 0.03, 0.025])


def test_carry_formula():
    assert carry_return(0.04, 0.037, 1.0, 0.02) == pytest.approx((0.04 - 0.02) * 1.0)


def test_flat_curve_zero_rolldown():
    assert rolldown_return(CF, FLAT, 1.0) == pytest.approx(0.0, abs=1e-9)


def test_upward_curve_positive_rolldown():
    assert rolldown_return(CF, UP, 1.0) > 0


def test_steeper_curve_more_rolldown():
    assert rolldown_return(CF, STEEP, 1.0) > rolldown_return(CF, UP, 1.0)


def test_inverted_curve_negative_rolldown():
    assert rolldown_return(CF, INV, 1.0) < 0


def test_total_is_carry_plus_rolldown():
    c = carry_return(0.04, 0.037, 1.0, 0.02)
    r = rolldown_return(CF, UP, 1.0)
    assert total_carry_rolldown(CF, UP, 0.04, 1.0, 0.02) == pytest.approx(c + r, abs=1e-9)


def test_negative_carry_when_financing_exceeds_coupon():
    assert carry_return(0.02, 0.03, 1.0, 0.05) < 0


def test_validation():
    with pytest.raises(ValueError):
        rolldown_return(CF, UP, 0)
