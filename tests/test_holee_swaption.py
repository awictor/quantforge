"""Ho-Lee bond option, coupon-bond option, and swaption (Jamshidian)."""

import pytest

from quantforge import (
    holee_swaption, holee_bond_option, holee_coupon_bond_option,
)
from quantforge.holee import holee_zero_coupon_bond


THETA, SIG, R0 = 0.005, 0.01, 0.03
T0 = 1.0
K = 0.03
PAY = [2.0, 3.0, 4.0, 5.0]


def test_single_cashflow_reduces_to_zcb_option():
    Kb = holee_zero_coupon_bond(0.03, 4.0, THETA, SIG)
    one = holee_coupon_bond_option(R0, T0, [(5.0, 1.0)], Kb, THETA, SIG, True)
    ref = holee_bond_option(R0, T0, 5.0, Kb, THETA, SIG, True)
    assert one == pytest.approx(ref, abs=1e-9)


def test_swaption_parity():
    payer = holee_swaption(R0, T0, PAY, K, THETA, SIG, payer=True)
    recv = holee_swaption(R0, T0, PAY, K, THETA, SIG, payer=False)
    prev = T0
    ann = 0.0
    for ti in PAY:
        ann += (ti - prev) * holee_zero_coupon_bond(R0, ti, THETA, SIG)
        prev = ti
    sr = ((holee_zero_coupon_bond(R0, T0, THETA, SIG)
           - holee_zero_coupon_bond(R0, PAY[-1], THETA, SIG)) / ann)
    assert (payer - recv) == pytest.approx(ann * (sr - K), abs=1e-8)


def test_bond_option_put_call_parity():
    Kb = 0.9
    c = holee_bond_option(R0, T0, 5.0, Kb, THETA, SIG, True)
    p = holee_bond_option(R0, T0, 5.0, Kb, THETA, SIG, False)
    P_bond = holee_zero_coupon_bond(R0, 5.0, THETA, SIG)
    P_opt = holee_zero_coupon_bond(R0, T0, THETA, SIG)
    assert (c - p) == pytest.approx(P_bond - Kb * P_opt, abs=1e-10)


def test_both_swaptions_positive():
    payer = holee_swaption(R0, T0, PAY, K, THETA, SIG, payer=True)
    recv = holee_swaption(R0, T0, PAY, K, THETA, SIG, payer=False)
    assert payer > 0.0 and recv > 0.0


def test_bad_pay_times_raise():
    with pytest.raises(ValueError):
        holee_swaption(R0, T0, [0.5, 2.0], K, THETA, SIG)
