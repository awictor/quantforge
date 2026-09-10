"""Barrier Monte Carlo vs the Reiner-Rubinstein closed form (with dividend q)."""

import pytest

from quantforge import (
    OptionType,
    barrier_mc,
    barrier_option,
    call_price,
)
from quantforge.exotics import Barrier


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.25
Q = 0.03
B = R - Q  # carry with a continuous dividend yield

CASES = [
    ("down-out", Barrier.DOWN_OUT, 90.0),
    ("down-in", Barrier.DOWN_IN, 90.0),
    ("up-out", Barrier.UP_OUT, 130.0),
    ("up-in", Barrier.UP_IN, 130.0),
]


@pytest.mark.slow
@pytest.mark.parametrize("name,kind,H", CASES)
def test_barrier_mc_matches_closed_form_with_dividend(name, kind, H):
    cf = barrier_option(S, K, H, T, R, SIGMA, OptionType.CALL, kind, b=B)
    mc = barrier_mc(S, K, H, T, R, SIGMA, OptionType.CALL, name, b=B,
                    n_steps=150, n_paths=120_000, seed=11)
    # Within ~3 standard errors of the continuous-monitoring closed form.
    assert abs(mc.price - cf) < 3.0 * mc.std_error + 0.05


@pytest.mark.slow
def test_brownian_bridge_removes_knockout_overpricing():
    # Naive discrete monitoring misses between-step crossings and over-prices a
    # knock-out; the bridge correction pulls it down toward the closed form.
    cf = barrier_option(S, K, 90.0, T, R, SIGMA, OptionType.CALL,
                        Barrier.DOWN_OUT, b=B)
    bridge = barrier_mc(S, K, 90.0, T, R, SIGMA, OptionType.CALL, "down-out",
                        b=B, n_steps=50, n_paths=80_000, seed=1)
    naive = barrier_mc(S, K, 90.0, T, R, SIGMA, OptionType.CALL, "down-out",
                       b=B, n_steps=50, n_paths=80_000, seed=1,
                       brownian_bridge=False)
    assert naive.price > bridge.price
    assert abs(bridge.price - cf) < abs(naive.price - cf)


@pytest.mark.slow
def test_in_out_parity():
    # Knock-in + knock-out = vanilla (no rebate), same paths.
    van = call_price(S, K, T, R, SIGMA, b=B)
    ki = barrier_mc(S, K, 90.0, T, R, SIGMA, OptionType.CALL, "down-in", b=B,
                    n_steps=100, n_paths=100_000, seed=3)
    ko = barrier_mc(S, K, 90.0, T, R, SIGMA, OptionType.CALL, "down-out", b=B,
                    n_steps=100, n_paths=100_000, seed=3)
    assert ki.price + ko.price == pytest.approx(van, abs=0.1)


def test_bad_barrier_raises():
    with pytest.raises(ValueError):
        barrier_mc(S, K, 90.0, T, R, SIGMA, barrier="sideways", n_paths=100)


def test_bad_H_raises():
    with pytest.raises(ValueError):
        barrier_mc(S, K, -1.0, T, R, SIGMA, barrier="down-out", n_paths=100)
