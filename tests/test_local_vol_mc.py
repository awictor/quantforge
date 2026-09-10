"""Tests for the local-volatility Monte Carlo pricer."""

import math

import pytest

from quantforge import local_vol_mc, call_price, put_price, OptionType


@pytest.mark.slow
def test_flat_local_vol_reproduces_black_scholes():
    r = local_vol_mc(100, 100, 1.0, 0.05, lambda S, tau: 0.2,
                     OptionType.CALL, n_steps=100, n_paths=100_000, seed=1)
    exact = call_price(100, 100, 1.0, 0.05, 0.2)
    # Within a few Euler-discretization / Monte Carlo standard errors.
    assert abs(r.price - exact) < 3 * r.std_error + 0.05


@pytest.mark.slow
def test_flat_local_vol_put():
    r = local_vol_mc(100, 105, 0.5, 0.03, lambda S, tau: 0.3,
                     OptionType.PUT, n_steps=100, n_paths=80_000, seed=2)
    exact = put_price(100, 105, 0.5, 0.03, 0.3)
    assert abs(r.price - exact) < 3 * r.std_error + 0.05


def test_skewed_local_vol_positive():
    # A downward-sloping local vol (leverage) gives a finite, positive price.
    r = local_vol_mc(100, 100, 1.0, 0.05, lambda S, tau: 0.2 * (100 / S) ** 0.5,
                     OptionType.CALL, n_steps=40, n_paths=8_000, seed=3)
    assert r.price > 0 and math.isfinite(r.price)


def test_time_dependent_local_vol_runs():
    # Vol rising with time should still price finitely.
    r = local_vol_mc(100, 100, 1.0, 0.05, lambda S, tau: 0.15 + 0.1 * tau,
                     OptionType.CALL, n_steps=40, n_paths=8_000, seed=4)
    assert r.price > 0


def test_reproducible_with_seed():
    fn = lambda S, tau: 0.2
    a = local_vol_mc(100, 100, 0.5, 0.05, fn, n_steps=50, n_paths=5_000, seed=9)
    b = local_vol_mc(100, 100, 0.5, 0.05, fn, n_steps=50, n_paths=5_000, seed=9)
    assert a.price == b.price


def test_rejects_bad_steps():
    with pytest.raises(ValueError):
        local_vol_mc(100, 100, 1.0, 0.05, lambda S, tau: 0.2, n_steps=0)
