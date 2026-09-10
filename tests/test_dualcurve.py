"""OIS/LIBOR dual-curve discounting and basis calibration."""

import math

import pytest

from quantforge import (
    DiscountCurve,
    dual_forward_rate,
    dual_par_swap_rate,
    dual_swap_value,
    dual_float_leg_value,
    dual_calibrate_basis,
)


OIS = DiscountCurve.from_zero_rates([1, 2, 3, 4, 5],
                                    [0.02, 0.022, 0.024, 0.025, 0.026])
PROJ = DiscountCurve.from_zero_rates([1, 2, 3, 4, 5],
                                     [0.025, 0.027, 0.029, 0.030, 0.031])


def test_single_curve_reduces_to_discount_curve_par():
    c = DiscountCurve.from_zero_rates([1, 2, 3, 5], [0.03, 0.032, 0.034, 0.036])
    pay = [1, 2, 3]
    assert dual_par_swap_rate(c, c, pay) == pytest.approx(c.par_swap_rate(pay),
                                                          abs=1e-12)


def test_forward_rate_flat_curve():
    c = DiscountCurve.from_zero_rates([1, 5, 10], [0.03, 0.03, 0.03])
    # Simply-compounded forward of a 3% continuous flat curve.
    f = dual_forward_rate(c, 1, 2)
    assert f == pytest.approx(math.exp(0.03) - 1.0, abs=1e-6)


def test_basis_recovery_exact():
    mats = [1, 2, 3, 4, 5]
    true_basis = 0.0015
    pars = [dual_par_swap_rate(OIS, PROJ, [k + 1 for k in range(m)], true_basis)
            for m in mats]
    b, rmse = dual_calibrate_basis(OIS, PROJ, mats, pars)
    assert b == pytest.approx(true_basis, abs=1e-6)
    assert rmse < 1e-6


def test_payer_swap_zero_at_par():
    pay = [1, 2, 3, 4, 5]
    par = dual_par_swap_rate(OIS, PROJ, pay)
    assert dual_swap_value(OIS, PROJ, pay, par, payer=True) == pytest.approx(
        0.0, abs=1e-12)


def test_payer_plus_receiver_is_zero():
    pay = [1, 2, 3, 4, 5]
    p = dual_swap_value(OIS, PROJ, pay, 0.03, payer=True)
    r = dual_swap_value(OIS, PROJ, pay, 0.03, payer=False)
    assert (p + r) == pytest.approx(0.0, abs=1e-12)


def test_basis_raises_float_leg():
    pay = [1, 2, 3, 4, 5]
    f0 = dual_float_leg_value(OIS, PROJ, pay, 0.0)
    fb = dual_float_leg_value(OIS, PROJ, pay, 0.002)
    assert fb > f0


def test_projection_above_ois_raises_par_rate():
    # A projection curve above OIS forecasts higher forwards -> higher par rate.
    pay = [1, 2, 3, 4, 5]
    high = dual_par_swap_rate(OIS, PROJ, pay)
    same = dual_par_swap_rate(OIS, OIS, pay)
    assert high > same
