"""Forward rate agreement valuation."""

import math

import pytest

from quantforge import fra_forward_rate as fwd, fra_value as fv


def _flat(r=0.03):
    return lambda t: math.exp(-r * t)


def test_zero_value_at_fair_rate():
    P = _flat()
    f = fwd(P, 1, 1.5)
    assert abs(fv(P, f, 1, 1.5)) < 1e-15


def test_payer_sign_versus_contract():
    P = _flat()
    assert fv(P, 0.02, 1, 1.5) > 0        # forward > contract
    assert fv(P, 0.05, 1, 1.5) < 0        # forward < contract


def test_receiver_is_negative_payer():
    P = _flat()
    assert abs(fv(P, 0.04, 1, 1.5, payer=True) + fv(P, 0.04, 1, 1.5, payer=False)) < 1e-15


def test_notional_scales():
    P = _flat()
    assert abs(fv(P, 0.02, 1, 1.5, notional=1e6) - 1e6 * fv(P, 0.02, 1, 1.5)) < 1e-6


def test_upward_curve_forward_rises():
    Pu = lambda t: math.exp(-(0.02 + 0.004 * t) * t) if t > 0 else 1.0
    assert fwd(Pu, 3, 3.5) > fwd(Pu, 1, 1.5)


def test_validation():
    with pytest.raises(ValueError):
        fv(_flat(), 0.03, 1.5, 1.0)       # t1 >= t2
    with pytest.raises(ValueError):
        fwd(_flat(), 2.0, 2.0)            # zero accrual
