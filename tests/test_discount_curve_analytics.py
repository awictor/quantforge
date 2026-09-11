"""DiscountCurve analytics: instantaneous forward, annuity, swap value / DV01."""

import math

import pytest

from quantforge import DiscountCurve


def _flat(r=0.03):
    return DiscountCurve.from_zero_rates([1, 2, 3, 4, 5], [r] * 5)


PAY = [1, 2, 3, 4, 5]


def test_instantaneous_forward_flat_equals_zero_rate():
    c = _flat(0.03)
    assert c.instantaneous_forward(2.0) == pytest.approx(0.03, abs=1e-6)
    assert c.instantaneous_forward(4.5) == pytest.approx(0.03, abs=1e-6)


def test_par_swap_has_zero_value():
    c = _flat(0.03)
    par = c.par_swap_rate(PAY)
    assert c.swap_value(PAY, par, payer=True) == pytest.approx(0.0, abs=1e-12)


def test_annuity_matches_manual_sum():
    c = _flat(0.03)
    manual = sum((PAY[i] - (0 if i == 0 else PAY[i - 1])) * c.df(PAY[i])
                 for i in range(len(PAY)))
    assert c.annuity(PAY) == pytest.approx(manual, abs=1e-12)


def test_payer_dv01_negative_receiver_positive():
    c = _flat(0.03)
    par = c.par_swap_rate(PAY)
    assert c.swap_dv01(PAY, par, payer=True) < 0.0
    assert c.swap_dv01(PAY, par, payer=False) > 0.0


def test_dv01_matches_pv01_scale():
    # For a par swap the DV01 magnitude is close to annuity * bump (the fixed leg
    # dominates the first-order change).
    c = _flat(0.03)
    par = c.par_swap_rate(PAY)
    dv = c.swap_dv01(PAY, par, payer=True, bump=1e-4)
    # Same order of magnitude as annuity * 1bp.
    assert abs(dv) == pytest.approx(c.annuity(PAY) * 1e-4, rel=0.5)


def test_bad_forward_time_raises():
    with pytest.raises(ValueError):
        _flat().instantaneous_forward(-1.0)
