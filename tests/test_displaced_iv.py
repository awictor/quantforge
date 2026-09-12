"""Displaced-diffusion implied volatility."""

import pytest

from quantforge import displaced_diffusion_price as dd, displaced_diffusion_implied_vol as ddiv
from quantforge.implied import implied_volatility


def test_round_trip_across_grid():
    S, K, t, r = 100, 100, 1.0, 0.05
    for sig in (0.1, 0.2, 0.35):
        for sh in (0.0, 25.0, 60.0):
            p = dd(S, K, t, r, sig, shift=sh)
            assert abs(ddiv(p, S, K, t, r, shift=sh) - sig) < 1e-6


def test_zero_shift_matches_bsm_iv():
    S, K, t, r = 100, 100, 1.0, 0.05
    p = dd(S, K, t, r, 0.25, shift=0.0)
    assert abs(ddiv(p, S, K, t, r, shift=0.0) - implied_volatility(p, S, K, t, r)) < 1e-6


def test_put_round_trip():
    p = dd(100, 90, 1.0, 0.05, 0.3, shift=40, option_type="put")
    assert abs(ddiv(p, 100, 90, 1.0, 0.05, shift=40, option_type="put") - 0.3) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        ddiv(1e6, 100, 100, 1.0, 0.05, shift=25)     # out of band
    with pytest.raises(ValueError):
        ddiv(5, 100, 100, 0.0, 0.05)                 # expiry
