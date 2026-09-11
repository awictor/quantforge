"""Variance-swap strike and VIX replicated from a SABR smile."""

import pytest

from quantforge import sabr_variance_swap_strike, sabr_vix, sabr_vol


F, T, R = 100.0, 1.0, 0.03


def test_flat_lognormal_variance():
    # beta = 1, nu -> 0: SABR reduces to a flat lognormal vol = alpha, so the
    # variance-swap strike is alpha^2.
    al = 0.2
    vs = sabr_variance_swap_strike(F, T, R, al, 1.0, 0.0, 1e-8,
                                   n_strikes=801, width=10.0)
    assert vs == pytest.approx(al * al, abs=2e-3)


def test_flat_vix_is_hundred_alpha():
    al = 0.2
    vx = sabr_vix(F, T, R, al, 1.0, 0.0, 1e-8, n_strikes=401, width=8.0)
    assert vx == pytest.approx(100.0 * al, abs=0.2)


def test_skew_variance_above_atm():
    vs = sabr_variance_swap_strike(F, T, R, 0.2, 0.5, -0.3, 0.4)
    atm = sabr_vol(F, F, T, 0.2, 0.5, -0.3, 0.4)
    assert vs > atm * atm


def test_positive():
    assert sabr_variance_swap_strike(F, T, R, 0.2, 0.5, -0.3, 0.4) > 0.0
    assert sabr_vix(F, T, R, 0.2, 0.5, -0.3, 0.4) > 0.0


def test_higher_volvol_raises_strike():
    lo = sabr_variance_swap_strike(F, T, R, 0.2, 0.5, -0.3, 0.2)
    hi = sabr_variance_swap_strike(F, T, R, 0.2, 0.5, -0.3, 0.6)
    assert hi > lo
