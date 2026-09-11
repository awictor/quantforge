"""European swaption under Vasicek via the coupon-bond-option identity."""

import pytest

from quantforge import vasicek_swaption
from quantforge.vasicek import zero_coupon_bond


KAPPA, THETA, SIG, R0 = 0.15, 0.05, 0.01, 0.04
T0 = 1.0
K = 0.03
PAY = [2.0, 3.0, 4.0, 5.0]


def _annuity_and_swap_rate():
    prev = T0
    ann = 0.0
    for ti in PAY:
        ann += (ti - prev) * zero_coupon_bond(R0, ti, KAPPA, THETA, SIG)
        prev = ti
    swap_rate = ((zero_coupon_bond(R0, T0, KAPPA, THETA, SIG)
                  - zero_coupon_bond(R0, PAY[-1], KAPPA, THETA, SIG)) / ann)
    return ann, swap_rate


def test_swaption_parity():
    payer = vasicek_swaption(R0, T0, PAY, K, KAPPA, THETA, SIG, payer=True)
    recv = vasicek_swaption(R0, T0, PAY, K, KAPPA, THETA, SIG, payer=False)
    ann, swap_rate = _annuity_and_swap_rate()
    assert (payer - recv) == pytest.approx(ann * (swap_rate - K), abs=1e-8)


def test_both_positive():
    payer = vasicek_swaption(R0, T0, PAY, K, KAPPA, THETA, SIG, payer=True)
    recv = vasicek_swaption(R0, T0, PAY, K, KAPPA, THETA, SIG, payer=False)
    assert payer > 0.0 and recv > 0.0


def test_atm_payer_equals_receiver():
    # At the forward swap rate the payer and receiver swaptions are worth the
    # same (parity's intrinsic term vanishes).
    _ann, swap_rate = _annuity_and_swap_rate()
    payer = vasicek_swaption(R0, T0, PAY, swap_rate, KAPPA, THETA, SIG, payer=True)
    recv = vasicek_swaption(R0, T0, PAY, swap_rate, KAPPA, THETA, SIG,
                            payer=False)
    assert payer == pytest.approx(recv, abs=1e-8)


def test_higher_strike_lowers_payer():
    lo = vasicek_swaption(R0, T0, PAY, 0.02, KAPPA, THETA, SIG, payer=True)
    hi = vasicek_swaption(R0, T0, PAY, 0.04, KAPPA, THETA, SIG, payer=True)
    assert hi < lo


def test_bad_pay_times_raise():
    with pytest.raises(ValueError):
        vasicek_swaption(R0, T0, [0.5, 2.0], K, KAPPA, THETA, SIG)
