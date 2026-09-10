"""Tests for Bachelier swaptions."""

import math

import pytest

from quantforge import (
    CapletPeriod, swaption_price, swaption_parity, annuity, bachelier_price,
    OptionType,
)


def _swap(n=5, rate=0.03):
    return [CapletPeriod(forward=rate, expiry=float(i), accrual=1.0,
                         discount=math.exp(-0.03 * i), sigma_n=0.01)
            for i in range(1, n + 1)]


def test_payer_receiver_parity():
    periods = _swap()
    sr, K = 0.032, 0.03
    payer = swaption_price(sr, K, 2.0, 0.01, periods, payer=True)
    receiver = swaption_price(sr, K, 2.0, 0.01, periods, payer=False)
    assert payer - receiver == pytest.approx(swaption_parity(sr, K, periods), abs=1e-12)


def test_value_is_annuity_times_bachelier():
    periods = _swap()
    sr, K, T, vol = 0.03, 0.035, 3.0, 0.012
    v = swaption_price(sr, K, T, vol, periods, payer=True)
    expected = annuity(periods) * bachelier_price(sr, K, T, 0.0, vol, OptionType.CALL)
    assert v == pytest.approx(expected, abs=1e-12)


def test_payer_increases_with_swap_rate():
    periods = _swap()
    lo = swaption_price(0.030, 0.035, 2.0, 0.01, periods, payer=True)
    hi = swaption_price(0.045, 0.035, 2.0, 0.01, periods, payer=True)
    assert hi > lo


def test_receiver_decreases_with_swap_rate():
    periods = _swap()
    lo = swaption_price(0.030, 0.035, 2.0, 0.01, periods, payer=False)
    hi = swaption_price(0.045, 0.035, 2.0, 0.01, periods, payer=False)
    assert lo > hi


def test_higher_vol_higher_value():
    periods = _swap()
    lo = swaption_price(0.03, 0.03, 2.0, 0.005, periods, payer=True)
    hi = swaption_price(0.03, 0.03, 2.0, 0.02, periods, payer=True)
    assert hi > lo


def test_handles_negative_swap_rate():
    periods = _swap(rate=-0.002)
    v = swaption_price(-0.002, 0.0, 1.0, 0.01, periods, payer=False)
    assert v > 0


def test_zero_expiry_is_annuity_intrinsic():
    periods = _swap()
    v = swaption_price(0.04, 0.03, 0.0, 0.01, periods, payer=True)
    assert v == pytest.approx(annuity(periods) * (0.04 - 0.03), abs=1e-12)


def test_annuity_is_sum_of_accrual_discount():
    periods = _swap()
    assert annuity(periods) == pytest.approx(sum(p.accrual * p.discount for p in periods))
