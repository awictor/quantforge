"""Corridor variance-swap fair strike replicated from a smile."""

import pytest

from quantforge import (
    corridor_variance_swap_from_smile as corr,
    variance_swap_from_smile,
)


S0, T, R, SIG = 100.0, 1.0, 0.03, 0.2


def _flat(K):
    return SIG


def test_narrower_corridor_accrues_less_variance():
    wide = corr(S0, T, R, _flat, 70, 140, n_strikes=400)
    narrow = corr(S0, T, R, _flat, 90, 110, n_strikes=400)
    assert narrow < wide


def test_nested_corridors_monotone_fixed_density():
    def kc(L, U):
        n = max(int(400 * (U - L) / 40), 50)   # keep dK roughly constant
        return corr(S0, T, R, _flat, L, U, n_strikes=n)
    vals = [kc(95, 105), kc(85, 120), kc(70, 140), kc(55, 170)]
    assert all(vals[i] <= vals[i + 1] + 1e-4 for i in range(len(vals) - 1))


def test_wide_corridor_near_full_variance_swap():
    full = variance_swap_from_smile(S0, T, R, _flat, n_strikes=401, width=8.0)
    wide = corr(S0, T, R, _flat, 40, 250, n_strikes=600)
    assert wide == pytest.approx(full, abs=5e-3)


def test_corridor_positive():
    assert corr(S0, T, R, _flat, 80, 130, n_strikes=300) > 0.0


def test_bad_corridor_raises():
    with pytest.raises(ValueError):
        corr(S0, T, R, _flat, 120, 80)      # lower > upper
    with pytest.raises(ValueError):
        corr(S0, 0.0, R, _flat, 80, 120)    # bad tenor
