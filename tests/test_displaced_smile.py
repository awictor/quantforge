"""BSM implied-vol smile produced by displaced diffusion (displaced_diffusion_smile)."""

import math

import pytest

from quantforge import (
    displaced_diffusion_smile, displaced_diffusion_price, implied_volatility,
)


S, T, R, SIG = 100.0, 1.0, 0.05, 0.2
STRIKES = [80.0, 90.0, 100.0, 110.0, 120.0]


def test_zero_shift_is_flat_at_sigma():
    sm = displaced_diffusion_smile(S, STRIKES, T, R, SIG, 0.0)
    assert all(v == pytest.approx(SIG, abs=1e-6) for _, v in sm)


def test_atm_vol_near_sigma():
    F = S * math.exp(R * T)
    sm = displaced_diffusion_smile(S, [F], T, R, SIG, 50.0)
    assert sm[0][1] == pytest.approx(SIG, abs=5e-3)


def test_positive_shift_downward_skew():
    sm = dict(displaced_diffusion_smile(S, [80.0, 120.0], T, R, SIG, 100.0))
    left = [v for k, v in sm.items() if k < 0][0]
    right = [v for k, v in sm.items() if k > 0][0]
    assert left > right


def test_larger_shift_steeper_skew():
    def slope(shift):
        sm = dict(displaced_diffusion_smile(S, [80.0, 120.0], T, R, SIG, shift))
        ks = sorted(sm)
        return (sm[ks[-1]] - sm[ks[0]]) / (ks[-1] - ks[0])
    assert slope(200.0) < slope(50.0) < slope(0.0) + 1e-9


def test_output_sorted_and_paired():
    sm = displaced_diffusion_smile(S, list(reversed(STRIKES)), T, R, SIG, 50.0)
    ks = [k for k, _ in sm]
    assert ks == sorted(ks)
    assert len(sm) == len(STRIKES)


def test_reprices_to_displaced():
    shift = 75.0
    F = S * math.exp(R * T)
    for lm, iv in displaced_diffusion_smile(S, STRIKES, T, R, SIG, shift):
        K = F * math.exp(lm)
        dd = displaced_diffusion_price(S, K, T, R, SIG, shift, "call")
        recov = implied_volatility(dd, S, K, T, R, "call", b=R)
        assert iv == pytest.approx(recov, abs=1e-9)
