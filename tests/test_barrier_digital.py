"""Tests for the barrier-contingent digital Monte Carlo pricer."""

import math

import pytest

from quantforge import barrier_digital_mc, cash_or_nothing, OptionType


@pytest.mark.slow
def test_knock_in_plus_knock_out_equals_plain_digital():
    # Touched-or-not partitions the paths, so KI + KO = the plain digital.
    S, K, H = 100, 100, 120
    ki = barrier_digital_mc(S, K, H, 1.0, 0.05, 0.2, OptionType.CALL, "up-in",
                            n_steps=200, n_paths=60_000, seed=1)
    ko = barrier_digital_mc(S, K, H, 1.0, 0.05, 0.2, OptionType.CALL, "up-out",
                            n_steps=200, n_paths=60_000, seed=1)
    plain = cash_or_nothing(S, K, 1.0, 0.05, 0.2, OptionType.CALL, cash=1.0)
    assert ki.price + ko.price == pytest.approx(plain, abs=3 * (ki.std_error + ko.std_error) + 1e-3)


def test_price_within_discounted_cash():
    r = 0.05
    res = barrier_digital_mc(100, 100, 120, 1.0, r, 0.2, OptionType.CALL, "up-in",
                             n_steps=100, n_paths=20_000, seed=2)
    assert 0.0 <= res.price <= math.exp(-r * 1.0) + 1e-9


@pytest.mark.slow
def test_far_knock_in_is_cheap():
    # A very high up-in barrier is rarely touched -> the KI digital is cheap.
    near = barrier_digital_mc(100, 100, 105, 1.0, 0.05, 0.3, OptionType.CALL,
                              "up-in", n_steps=100, n_paths=20_000, seed=3)
    far = barrier_digital_mc(100, 100, 160, 1.0, 0.05, 0.3, OptionType.CALL,
                             "up-in", n_steps=100, n_paths=20_000, seed=3)
    assert far.price < near.price


def test_down_in_put():
    res = barrier_digital_mc(100, 100, 80, 1.0, 0.05, 0.3, OptionType.PUT,
                             "down-in", n_steps=100, n_paths=20_000, seed=4)
    assert res.price > 0


def test_cash_scales_linearly():
    a = barrier_digital_mc(100, 100, 120, 1.0, 0.05, 0.2, OptionType.CALL,
                           "up-in", cash=1.0, n_steps=50, n_paths=10_000, seed=5)
    b = barrier_digital_mc(100, 100, 120, 1.0, 0.05, 0.2, OptionType.CALL,
                           "up-in", cash=10.0, n_steps=50, n_paths=10_000, seed=5)
    assert b.price == pytest.approx(10.0 * a.price, rel=1e-9)


def test_reproducible_with_seed():
    a = barrier_digital_mc(100, 100, 120, 1.0, 0.05, 0.2, seed=9,
                           n_steps=50, n_paths=5_000)
    b = barrier_digital_mc(100, 100, 120, 1.0, 0.05, 0.2, seed=9,
                           n_steps=50, n_paths=5_000)
    assert a.price == b.price


def test_rejects_bad_barrier():
    with pytest.raises(ValueError):
        barrier_digital_mc(100, 100, 120, 1.0, 0.05, 0.2, barrier="sideways",
                           n_paths=100)
