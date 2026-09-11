"""Range binary (double digital) and supershare options (exotics)."""

import math

import pytest

from quantforge import (
    range_binary, supershare, cash_or_nothing, asset_or_nothing,
)


S, R, SIG, T = 100.0, 0.05, 0.2, 1.0
KL, KH = 95.0, 110.0


def test_range_binary_equals_cash_or_nothing_difference():
    rb = range_binary(S, KL, KH, T, R, SIG)
    diff = (cash_or_nothing(S, KL, T, R, SIG, "call")
            - cash_or_nothing(S, KH, T, R, SIG, "call"))
    assert rb == pytest.approx(diff, abs=1e-12)


def test_supershare_equals_asset_or_nothing_difference():
    ss = supershare(S, KL, KH, T, R, SIG)
    diff = (asset_or_nothing(S, KL, T, R, SIG, "call")
            - asset_or_nothing(S, KH, T, R, SIG, "call")) / KL
    assert ss == pytest.approx(diff, abs=1e-12)


def test_range_binary_matches_monte_carlo():
    import random
    rng = random.Random(3)
    N = 500000
    disc = math.exp(-R * T)
    hits = 0
    for _ in range(N):
        z = rng.gauss(0.0, 1.0)
        ST = S * math.exp((R - 0.5 * SIG * SIG) * T + SIG * math.sqrt(T) * z)
        if KL <= ST <= KH:
            hits += 1
    mc = disc * hits / N
    assert range_binary(S, KL, KH, T, R, SIG) == pytest.approx(mc, rel=0.02)


def test_supershare_matches_monte_carlo():
    import random
    rng = random.Random(3)
    N = 500000
    disc = math.exp(-R * T)
    acc = 0.0
    for _ in range(N):
        z = rng.gauss(0.0, 1.0)
        ST = S * math.exp((R - 0.5 * SIG * SIG) * T + SIG * math.sqrt(T) * z)
        if KL <= ST <= KH:
            acc += ST / KL
    mc = disc * acc / N
    assert supershare(S, KL, KH, T, R, SIG) == pytest.approx(mc, rel=0.02)


def test_range_binary_bounded_by_cash():
    # Probability-weighted cash payoff can never exceed the discounted cash.
    rb = range_binary(S, KL, KH, T, R, SIG, cash=1.0)
    assert 0.0 < rb < math.exp(-R * T)


def test_wide_corridor_approaches_full_cash():
    # A very wide corridor collects essentially the whole discounted cash.
    rb = range_binary(S, 1e-6, 1e6, T, R, SIG, cash=1.0)
    assert rb == pytest.approx(math.exp(-R * T), rel=1e-6)


def test_zero_vol_pays_when_forward_inside():
    fwd = S * math.exp(R * T)
    inside = range_binary(S, fwd - 1.0, fwd + 1.0, T, R, 0.0)
    outside = range_binary(S, fwd + 5.0, fwd + 10.0, T, R, 0.0)
    assert inside == pytest.approx(math.exp(-R * T), abs=1e-12)
    assert outside == pytest.approx(0.0, abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        range_binary(S, 110.0, 95.0, T, R, SIG)  # K_high <= K_low
    with pytest.raises(ValueError):
        supershare(S, 110.0, 95.0, T, R, SIG)
