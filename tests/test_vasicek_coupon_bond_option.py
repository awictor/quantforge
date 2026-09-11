"""European coupon-bond option under Vasicek via Jamshidian (vasicek_coupon_bond_option)."""

import pytest

from quantforge import vasicek_coupon_bond_option, OptionType
from quantforge.vasicek import bond_option, zero_coupon_bond


KAPPA, THETA, SIG, R0 = 0.15, 0.05, 0.01, 0.04
TOPT = 1.0
CFS = [(2.0, 0.03), (3.0, 0.03), (4.0, 0.03), (5.0, 1.03)]
KC = 0.95


def test_single_cashflow_reduces_to_zcb_option():
    K = zero_coupon_bond(0.05, 4.0, KAPPA, THETA, SIG)
    one = vasicek_coupon_bond_option(R0, TOPT, [(5.0, 1.0)], K, KAPPA, THETA,
                                     SIG, OptionType.CALL)
    ref = bond_option(R0, TOPT, 5.0, K, KAPPA, THETA, SIG, OptionType.CALL)
    assert one == pytest.approx(ref, abs=1e-9)


def test_put_call_parity():
    c = vasicek_coupon_bond_option(R0, TOPT, CFS, KC, KAPPA, THETA, SIG,
                                   OptionType.CALL)
    p = vasicek_coupon_bond_option(R0, TOPT, CFS, KC, KAPPA, THETA, SIG,
                                   OptionType.PUT)
    fwd_less_strike = (sum(c_ * zero_coupon_bond(R0, ti, KAPPA, THETA, SIG)
                           for ti, c_ in CFS)
                       - KC * zero_coupon_bond(R0, TOPT, KAPPA, THETA, SIG))
    assert (c - p) == pytest.approx(fwd_less_strike, abs=1e-8)


def test_prices_positive():
    c = vasicek_coupon_bond_option(R0, TOPT, CFS, KC, KAPPA, THETA, SIG,
                                   OptionType.CALL)
    p = vasicek_coupon_bond_option(R0, TOPT, CFS, KC, KAPPA, THETA, SIG,
                                   OptionType.PUT)
    assert c > 0.0 and p > 0.0


def test_higher_strike_lowers_call():
    lo = vasicek_coupon_bond_option(R0, TOPT, CFS, 0.93, KAPPA, THETA, SIG,
                                    OptionType.CALL)
    hi = vasicek_coupon_bond_option(R0, TOPT, CFS, 0.97, KAPPA, THETA, SIG,
                                    OptionType.CALL)
    assert hi < lo


def test_bad_cashflow_time_raises():
    with pytest.raises(ValueError):
        vasicek_coupon_bond_option(R0, TOPT, [(0.5, 1.0)], KC, KAPPA, THETA, SIG)
