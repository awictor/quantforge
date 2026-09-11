"""CIR coupon-bond option and swaption via Jamshidian (cir_swaption)."""

import pytest

from quantforge import (
    cir_swaption, cir_coupon_bond_option, cir_bond_option, OptionType,
)
from quantforge.cir import cir_zero_coupon_bond


KAPPA, THETA, SIG, R0 = 0.3, 0.04, 0.05, 0.04
T0 = 1.0
K = 0.04
PAY = [2.0, 3.0, 4.0, 5.0]


def test_single_cashflow_reduces_to_zcb_option():
    Kb = cir_zero_coupon_bond(0.04, 4.0, KAPPA, THETA, SIG)
    one = cir_coupon_bond_option(R0, T0, [(5.0, 1.0)], Kb, KAPPA, THETA, SIG,
                                 OptionType.CALL)
    ref = cir_bond_option(R0, T0, 5.0, Kb, KAPPA, THETA, SIG, OptionType.CALL)
    assert one == pytest.approx(ref, abs=1e-9)


def test_swaption_parity():
    payer = cir_swaption(R0, T0, PAY, K, KAPPA, THETA, SIG, payer=True)
    recv = cir_swaption(R0, T0, PAY, K, KAPPA, THETA, SIG, payer=False)
    prev = T0
    ann = 0.0
    for ti in PAY:
        ann += (ti - prev) * cir_zero_coupon_bond(R0, ti, KAPPA, THETA, SIG)
        prev = ti
    sr = ((cir_zero_coupon_bond(R0, T0, KAPPA, THETA, SIG)
           - cir_zero_coupon_bond(R0, PAY[-1], KAPPA, THETA, SIG)) / ann)
    assert (payer - recv) == pytest.approx(ann * (sr - K), abs=1e-7)


def test_both_positive():
    payer = cir_swaption(R0, T0, PAY, K, KAPPA, THETA, SIG, payer=True)
    recv = cir_swaption(R0, T0, PAY, K, KAPPA, THETA, SIG, payer=False)
    assert payer > 0.0 and recv > 0.0


def test_higher_strike_lowers_payer():
    lo = cir_swaption(R0, T0, PAY, 0.03, KAPPA, THETA, SIG, payer=True)
    hi = cir_swaption(R0, T0, PAY, 0.05, KAPPA, THETA, SIG, payer=True)
    assert hi < lo


def test_bad_pay_times_raise():
    with pytest.raises(ValueError):
        cir_swaption(R0, T0, [0.5, 2.0], K, KAPPA, THETA, SIG)
