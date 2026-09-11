"""Ikeda-Kunitomo double-barrier knock-out call (double_knock_out_call)."""

import math

import pytest

from quantforge import double_knock_out_call
from quantforge.bsm import call_price


S, K, T, R, SIG = 100.0, 95.0, 0.5, 0.05, 0.25


def test_wide_barriers_recover_vanilla():
    # Barriers far from spot are almost never touched: price -> vanilla call.
    v = double_knock_out_call(S, K, 40.0, 400.0, T, R, SIG)
    assert v == pytest.approx(call_price(S, K, T, R, SIG), rel=1e-6)


def test_tighter_band_lowers_value():
    wide = double_knock_out_call(S, K, 80.0, 130.0, T, R, SIG)
    narrow = double_knock_out_call(S, K, 90.0, 112.0, T, R, SIG)
    assert 0.0 <= narrow < wide


def test_higher_vol_lowers_value():
    # More vol -> more likely to hit a barrier -> lower knock-out value.
    lo = double_knock_out_call(S, K, 88.0, 116.0, T, R, 0.15)
    hi = double_knock_out_call(S, K, 88.0, 116.0, T, R, 0.35)
    assert hi < lo


def test_price_nonnegative_and_below_vanilla():
    v = double_knock_out_call(S, K, 88.0, 116.0, T, R, SIG)
    assert 0.0 <= v <= call_price(S, K, T, R, SIG)


@pytest.mark.slow
def test_matches_fine_step_monte_carlo():
    import random

    L, U, sig = 90.0, 115.0, 0.25
    closed = double_knock_out_call(S, K, L, U, T, R, sig)
    rng = random.Random(7)
    steps, paths = 6000, 40000
    dt = T / steps
    drift = (R - 0.5 * sig * sig) * dt
    vol = sig * math.sqrt(dt)
    disc = math.exp(-R * T)
    acc = 0.0
    for _ in range(paths):
        x = S
        alive = True
        for _ in range(steps):
            x *= math.exp(drift + vol * rng.gauss(0.0, 1.0))
            if x <= L or x >= U:
                alive = False
                break
        if alive:
            acc += max(x - K, 0.0)
    mc = disc * acc / paths
    # Discrete monitoring overstates survival (bias decays like sqrt(dt)), so the
    # MC estimate sits slightly above the continuously-monitored closed form.
    assert closed <= mc <= closed + 0.1


def test_requires_ordering():
    with pytest.raises(ValueError):
        double_knock_out_call(S, K, 120.0, 90.0, T, R, SIG)
    with pytest.raises(ValueError):
        double_knock_out_call(S, K, 100.0, 120.0, T, R, SIG)  # S == L
