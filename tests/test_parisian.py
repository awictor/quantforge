"""Tests for the Parisian barrier Monte Carlo pricer."""

import pytest

from quantforge import (
    parisian_barrier_mc, barrier_option, call_price, Barrier, OptionType,
)


@pytest.mark.slow
def test_knock_in_plus_knock_out_is_vanilla():
    S, K, H = 100, 100, 90
    ko = parisian_barrier_mc(S, K, H, 1.0, 0.05, 0.25, window=0.1,
                             option_type=OptionType.CALL, barrier="down-out",
                             n_steps=252, n_paths=40_000, seed=1)
    ki = parisian_barrier_mc(S, K, H, 1.0, 0.05, 0.25, window=0.1,
                             option_type=OptionType.CALL, barrier="down-in",
                             n_steps=252, n_paths=40_000, seed=1)
    van = call_price(S, K, 1.0, 0.05, 0.25)
    assert ko.price + ki.price == pytest.approx(van, abs=3 * (ko.std_error + ki.std_error) + 1e-3)


@pytest.mark.slow
def test_parisian_knockout_above_standard_knockout():
    # A Parisian barrier is harder to activate than an instantaneous one, so a
    # Parisian knock-out (holder keeps the option unless activated) is worth
    # more than the standard knock-out.
    S, K, H = 100, 100, 90
    par = parisian_barrier_mc(S, K, H, 1.0, 0.05, 0.25, window=0.1,
                              option_type=OptionType.CALL, barrier="down-out",
                              n_steps=252, n_paths=40_000, seed=2)
    std = barrier_option(S, K, H, 1.0, 0.05, 0.25, OptionType.CALL, Barrier.DOWN_OUT)
    assert par.price > std


@pytest.mark.slow
def test_longer_window_raises_knockout_value():
    S, K, H = 100, 100, 90
    short = parisian_barrier_mc(S, K, H, 1.0, 0.05, 0.25, window=0.1,
                                barrier="down-out", n_steps=252, n_paths=40_000, seed=3)
    long_w = parisian_barrier_mc(S, K, H, 1.0, 0.05, 0.25, window=0.3,
                                 barrier="down-out", n_steps=252, n_paths=40_000, seed=3)
    assert long_w.price > short.price


def test_reproducible_with_seed():
    a = parisian_barrier_mc(100, 100, 90, 0.5, 0.05, 0.25, window=0.1,
                            n_steps=50, n_paths=4_000, seed=9)
    b = parisian_barrier_mc(100, 100, 90, 0.5, 0.05, 0.25, window=0.1,
                            n_steps=50, n_paths=4_000, seed=9)
    assert a.price == b.price


def test_rejects_bad_window():
    with pytest.raises(ValueError):
        parisian_barrier_mc(100, 100, 90, 1.0, 0.05, 0.25, window=0.0, n_paths=100)
    with pytest.raises(ValueError):
        parisian_barrier_mc(100, 100, 90, 1.0, 0.05, 0.25, window=2.0, n_paths=100)


def test_rejects_bad_barrier():
    with pytest.raises(ValueError):
        parisian_barrier_mc(100, 100, 90, 1.0, 0.05, 0.25, window=0.1,
                            barrier="sideways", n_paths=100)
