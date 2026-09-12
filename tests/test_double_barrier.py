"""Double-barrier knock-out call (Kunitomo-Ikeda closed form)."""

import math

import pytest

from quantforge import double_knockout_call, call_price

S, K, T, R, SIG = 100.0, 100.0, 1.0, 0.05, 0.25
B = R


def test_below_vanilla():
    van = call_price(S, K, T, R, SIG, b=B)
    assert 0.0 <= double_knockout_call(S, K, 80, 130, T, R, SIG, b=B) <= van


def test_far_barriers_approach_vanilla():
    van = call_price(S, K, T, R, SIG, b=B)
    far = double_knockout_call(S, K, 1.0, 1e5, T, R, SIG, b=B)
    assert abs(far - van) < 0.05


def test_series_converges_in_terms():
    lo = double_knockout_call(S, K, 80, 130, T, R, SIG, b=B, terms=4)
    hi = double_knockout_call(S, K, 80, 130, T, R, SIG, b=B, terms=12)
    assert abs(lo - hi) < 1e-6


def test_tighter_corridor_lowers_value():
    tight = double_knockout_call(S, K, 90, 115, T, R, SIG, b=B)
    wide = double_knockout_call(S, K, 70, 140, T, R, SIG, b=B)
    assert tight < wide


@pytest.mark.slow
def test_matches_monte_carlo_in_continuous_limit():
    # Discrete MC monitoring overstates the price; as the step count rises it
    # converges down toward the continuous-monitoring closed form.
    import random

    cf = double_knockout_call(S, K, 80, 130, T, R, SIG, b=B)
    L, U, M, N = 80.0, 130.0, 6000, 30000
    random.seed(11)
    dt = T / M
    drift = (B - 0.5 * SIG * SIG) * dt
    vol = SIG * math.sqrt(dt)
    disc = math.exp(-R * T)
    acc = 0.0
    for _ in range(N):
        s = S
        alive = True
        for _ in range(M):
            s *= math.exp(drift + vol * random.gauss(0, 1))
            if s <= L or s >= U:
                alive = False
                break
        if alive:
            acc += max(s - K, 0.0)
    mc = disc * acc / N
    # MC (discrete) should sit above the closed form but within ~0.15 at M=6000.
    assert cf < mc
    assert abs(mc - cf) < 0.15


def test_validation():
    with pytest.raises(ValueError):
        double_knockout_call(S, K, 110, 130, T, R, SIG)  # L >= S
    with pytest.raises(ValueError):
        double_knockout_call(S, K, 80, 90, T, R, SIG)    # U <= S
    with pytest.raises(ValueError):
        double_knockout_call(S, K, 80, 130, 0, R, SIG)   # t = 0
