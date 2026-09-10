"""Tests for the average-strike Asian Monte Carlo pricer."""

import pytest

from quantforge import average_strike_asian_mc, OptionType


def test_call_and_put_positive():
    c = average_strike_asian_mc(100, 1.0, 0.05, 0.3, OptionType.CALL,
                                n_steps=50, n_paths=60_000, seed=1)
    p = average_strike_asian_mc(100, 1.0, 0.05, 0.3, OptionType.PUT,
                                n_steps=50, n_paths=60_000, seed=1)
    assert c.price > 0 and p.price > 0


def test_single_step_is_zero():
    # With one monitoring date the average equals the terminal spot, so the
    # payoff max(S_T - A, 0) is identically zero.
    z = average_strike_asian_mc(100, 1.0, 0.05, 0.3, OptionType.CALL,
                                n_steps=1, n_paths=10_000, seed=2)
    assert z.price == pytest.approx(0.0, abs=1e-9)


def test_higher_vol_raises_price():
    lo = average_strike_asian_mc(100, 1.0, 0.05, 0.3, n_steps=50, n_paths=40_000, seed=3)
    hi = average_strike_asian_mc(100, 1.0, 0.05, 0.6, n_steps=50, n_paths=40_000, seed=3)
    assert hi.price > lo.price


def test_reproducible_with_seed():
    a = average_strike_asian_mc(100, 0.5, 0.05, 0.25, n_steps=20, n_paths=5_000, seed=9)
    b = average_strike_asian_mc(100, 0.5, 0.05, 0.25, n_steps=20, n_paths=5_000, seed=9)
    assert a.price == b.price


def test_confidence_interval_brackets_price():
    r = average_strike_asian_mc(100, 1.0, 0.05, 0.3, n_steps=50, n_paths=20_000, seed=4)
    lo, hi = r.confidence_interval()
    assert lo < r.price < hi


def test_rejects_bad_steps():
    with pytest.raises(ValueError):
        average_strike_asian_mc(100, 1.0, 0.05, 0.3, n_steps=0)
