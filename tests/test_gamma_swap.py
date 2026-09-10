"""Gamma (price-weighted variance) swap strike from a smile."""

import math

import pytest

from quantforge import gamma_swap_from_smile, variance_swap_from_smile


S0, T, R = 100.0, 1.0, 0.03


def _flat(K):
    return 0.2


def test_flat_smile_gives_variance():
    g = gamma_swap_from_smile(S0, T, R, _flat, n_strikes=601, width=10.0)
    assert g == pytest.approx(0.04, abs=5e-3)


def test_denser_strip_improves_flat_variance():
    # At a fixed width, more strikes (finer dK) reduce the trapezoidal error.
    e_coarse = abs(gamma_swap_from_smile(S0, T, R, _flat, n_strikes=101,
                                         width=8.0) - 0.04)
    e_fine = abs(gamma_swap_from_smile(S0, T, R, _flat, n_strikes=801,
                                       width=8.0) - 0.04)
    assert e_fine < e_coarse


def test_downward_skew_gamma_below_variance():
    # S-weighting down-weights the low-strike high-vol puts, so a downward-skewed
    # gamma swap is worth less than the same smile's variance swap.
    def smile(K):
        return max(0.05, 0.2 + 0.15 * math.log(S0 / K))
    g = gamma_swap_from_smile(S0, T, R, smile, n_strikes=401, width=6.0)
    v = variance_swap_from_smile(S0, T, R, smile, n_strikes=401, width=6.0)
    assert g < v


def test_positive():
    assert gamma_swap_from_smile(S0, T, R, _flat, n_strikes=301, width=6.0) > 0.0


def test_bad_tenor_raises():
    with pytest.raises(ValueError):
        gamma_swap_from_smile(S0, 0.0, R, _flat)
