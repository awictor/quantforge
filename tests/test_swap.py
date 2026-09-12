"""Single-curve vanilla interest-rate swap valuation."""

import math

import pytest

from quantforge import (
    vanilla_swap_value as sv, single_curve_par_swap_rate as psr, swap_annuity as ann,
)


def _flat(r=0.03):
    return lambda t: math.exp(-r * t)


def test_par_rate_gives_zero_value():
    P = _flat()
    par = psr(P, 0, 5, freq=2)
    assert abs(sv(P, par, 0, 5, freq=2)) < 1e-12


def test_receiver_is_negative_payer():
    P = _flat()
    assert abs(sv(P, 0.04, 0, 5, payer=True) + sv(P, 0.04, 0, 5, payer=False)) < 1e-15


def test_lower_fixed_raises_payer_value():
    P = _flat()
    assert sv(P, 0.02, 0, 5) > sv(P, 0.05, 0, 5)


def test_notional_scales_linearly():
    P = _flat()
    assert abs(sv(P, 0.04, 0, 5, notional=1e6) - 1e6 * sv(P, 0.04, 0, 5)) < 1e-6


def test_forward_starting_swap_zero_at_par():
    P = _flat()
    par = psr(P, 1, 6, freq=2)
    assert abs(sv(P, par, 1, 6, freq=2)) < 1e-12


def test_annuity_positive():
    assert ann(_flat(), 0, 5) > 0


def test_validation():
    with pytest.raises(ValueError):
        psr(_flat(), 0, 0.1, freq=2)      # no full period
