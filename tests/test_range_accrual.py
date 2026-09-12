"""Range-accrual note (closed form)."""

import math

import pytest

from quantforge import range_accrual_note

S, T, R, SIG, CPN, M = 100.0, 1.0, 0.05, 0.25, 0.06, 12
B = R


def test_full_range_equals_discounted_coupon():
    full = range_accrual_note(S, 1e-6, 1e9, T, R, SIG, CPN, M, b=B)
    assert abs(full - CPN * math.exp(-R * T)) < 1e-6


def test_wider_band_raises_value():
    wide = range_accrual_note(S, 80, 125, T, R, SIG, CPN, M, b=B)
    narrow = range_accrual_note(S, 95, 106, T, R, SIG, CPN, M, b=B)
    assert wide > narrow


def test_higher_vol_lowers_value():
    lo_vol = range_accrual_note(S, 90, 110, T, R, SIG, CPN, M, b=B)
    hi_vol = range_accrual_note(S, 90, 110, T, R, 0.5, CPN, M, b=B)
    assert hi_vol < lo_vol


def test_nonnegative_and_bounded():
    pv = range_accrual_note(S, 90, 110, T, R, SIG, CPN, M, b=B)
    assert 0.0 <= pv <= CPN * math.exp(-R * T)


def test_notional_scales_linearly():
    one = range_accrual_note(S, 90, 110, T, R, SIG, CPN, M, b=B, notional=1.0)
    hundred = range_accrual_note(S, 90, 110, T, R, SIG, CPN, M, b=B, notional=100.0)
    assert abs(hundred - 100.0 * one) < 1e-12


@pytest.mark.slow
def test_matches_monte_carlo():
    import random

    pv = range_accrual_note(S, 90, 110, T, R, SIG, CPN, M, b=B)
    L, U, N = 90.0, 110.0, 200000
    random.seed(5)
    disc = math.exp(-R * T)
    acc = 0.0
    for _ in range(N):
        cnt = 0
        for i in range(1, M + 1):
            ti = T * i / M
            z = random.gauss(0, 1)
            sti = S * math.exp((B - 0.5 * SIG * SIG) * ti + SIG * math.sqrt(ti) * z)
            if L <= sti <= U:
                cnt += 1
        acc += CPN * (cnt / M)
    mc = disc * acc / N
    assert abs(mc - pv) < 5e-4


def test_validation():
    with pytest.raises(ValueError):
        range_accrual_note(S, 110, 90, T, R, SIG, CPN, M)   # L >= U
    with pytest.raises(ValueError):
        range_accrual_note(S, 90, 110, T, R, SIG, CPN, 0)   # observations < 1
    with pytest.raises(ValueError):
        range_accrual_note(S, 90, 110, 0, R, SIG, CPN, M)   # t = 0
