"""Digital (binary) degree-day option."""

import math

import pytest

from quantforge import degree_day_digital as dd


MU, K, SIG, R, T, PAY = 1200, 1000, 150, 0.04, 0.5, 100000
DISC = math.exp(-R * T)


def test_call_put_sum_to_discounted_payout():
    c = dd(MU, K, SIG, R, T, PAY, is_call=True)
    p = dd(MU, K, SIG, R, T, PAY, is_call=False)
    assert abs((c + p) - DISC * PAY) < 1e-6


def test_atm_is_half():
    assert abs(dd(MU, MU, SIG, R, T, PAY) - DISC * PAY * 0.5) < 1e-6


def test_deep_itm_call_full_payout():
    assert abs(dd(MU, 200, SIG, R, T, PAY) - DISC * PAY) < 1e-3


def test_monotone_in_strike():
    assert dd(MU, 900, SIG, R, T, PAY) > dd(MU, 1100, SIG, R, T, PAY)


def test_zero_vol_step():
    assert dd(1200, 1000, 0.0, R, T, PAY) == DISC * PAY
    assert dd(900, 1000, 0.0, R, T, PAY) == 0.0


def test_validation():
    with pytest.raises(ValueError):
        dd(MU, K, -1.0, R, T, PAY)
    with pytest.raises(ValueError):
        dd(MU, K, SIG, R, T, -1.0)
