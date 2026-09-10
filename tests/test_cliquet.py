"""Tests for the capped-cliquet Monte Carlo pricer."""

import math

import pytest

from quantforge import capped_cliquet_mc


RESETS = [0.25, 0.5, 0.75, 1.0]


def test_global_cap_bounds_price():
    # The discounted payoff can never exceed the discounted global cap.
    res = capped_cliquet_mc(100, 1.0, 0.05, 0.3, RESETS, local_cap=0.05,
                            global_cap=0.15, n_paths=20_000, seed=1)
    assert res.price <= 0.15 * math.exp(-0.05) + 3 * res.std_error


def test_price_nonnegative_with_zero_floors():
    res = capped_cliquet_mc(100, 1.0, 0.03, 0.25, RESETS, local_cap=0.04,
                            n_paths=20_000, seed=2)
    assert res.price >= 0.0


def test_tighter_local_cap_lowers_price():
    loose = capped_cliquet_mc(100, 1.0, 0.03, 0.3, RESETS, local_cap=0.10,
                              n_paths=20_000, seed=3)
    tight = capped_cliquet_mc(100, 1.0, 0.03, 0.3, RESETS, local_cap=0.02,
                              n_paths=20_000, seed=3)
    assert tight.price < loose.price


def test_global_cap_lowers_price():
    uncapped = capped_cliquet_mc(100, 1.0, 0.03, 0.3, RESETS, local_cap=0.05,
                                 n_paths=20_000, seed=4)
    capped = capped_cliquet_mc(100, 1.0, 0.03, 0.3, RESETS, local_cap=0.05,
                               global_cap=0.08, n_paths=20_000, seed=4)
    assert capped.price <= uncapped.price + 1e-9


def test_reproducible_with_seed():
    a = capped_cliquet_mc(100, 1.0, 0.03, 0.25, RESETS, local_cap=0.05,
                          n_paths=5_000, seed=42)
    b = capped_cliquet_mc(100, 1.0, 0.03, 0.25, RESETS, local_cap=0.05,
                          n_paths=5_000, seed=42)
    assert a.price == b.price
    assert a.std_error == b.std_error


def test_zero_vol_is_deterministic_capped_drift():
    # With zero vol each period return is exp(r*dt)-1; sum them, clip, discount.
    r, sigma = 0.05, 0.0
    res = capped_cliquet_mc(100, 1.0, r, sigma, RESETS, local_cap=None,
                            global_cap=None, n_paths=2, seed=1)
    dt = 0.25
    per = math.exp(r * dt) - 1.0
    expected = math.exp(-r * 1.0) * 4 * per
    assert res.price == pytest.approx(expected, abs=1e-9)


def test_rejects_bad_schedule():
    with pytest.raises(ValueError):
        capped_cliquet_mc(100, 1.0, 0.03, 0.25, [0.5, 0.5], n_paths=100, seed=1)
