"""Unified exponential-Levy implied-vol surface builder."""

import math

import pytest

from quantforge import LevySurface, levy_psi
from quantforge.nig import nig_smile
from quantforge.variancegamma import variance_gamma_smile


S, R, Q = 100.0, 0.03, 0.0


def test_surface_reprices_nig_smile():
    params = (18.0, -6.0, 0.55)
    surf = LevySurface("nig", params, S, R, Q)
    t = 0.5
    F = S * math.exp(R * t)
    for K, iv_ref in [(K, iv) for K, (_, iv) in
                      zip([80, 90, 100, 110, 120],
                          nig_smile(S, [80, 90, 100, 110, 120], t, R, *params))]:
        iv_s = surf.implied_vol(math.log(K / F), t)
        assert iv_s == pytest.approx(iv_ref, abs=1e-3)


def test_surface_reprices_vg_smile():
    params = (0.2, 0.35, -0.25)
    surf = LevySurface("vg", params, S, R, Q)
    t = 1.0
    F = S * math.exp(R * t)
    ref = variance_gamma_smile(S, [90, 100, 110], t, R, *params)
    for (k_ref, iv_ref), K in zip(ref, [90, 100, 110]):
        iv_s = surf.implied_vol(math.log(K / F), t)
        assert iv_s == pytest.approx(iv_ref, abs=1e-3)


def test_levy_surface_is_calendar_arbitrage_free():
    # A genuine Levy law has total variance non-decreasing in maturity.
    surf = LevySurface("nig", (18.0, -6.0, 0.55), S, R, Q)
    expiries = [0.1, 0.25, 0.5, 1.0, 2.0]
    assert surf.is_calendar_arbitrage_free(expiries)
    assert surf.calendar_violations(expiries) == []


def test_total_variance_increases_with_maturity():
    surf = LevySurface("meixner", (0.3, -0.3, 0.5), S, R, Q)
    w_short = surf.total_variance(0.0, 0.25)
    w_long = surf.total_variance(0.0, 1.0)
    assert w_long > w_short


def test_grid_shape_and_skew():
    surf = LevySurface("nig", (18.0, -6.0, 0.55), S, R, Q)
    g = surf.grid([0.25, 1.0], [-0.1, 0.0, 0.1])
    assert sorted(g) == [0.25, 1.0]
    assert len(g[0.25]) == 3
    # beta < 0 => downward skew at each expiry.
    for t in g:
        vols = [iv for _, iv in g[t]]
        assert vols[0] > vols[-1]


def test_levy_psi_factory_matches_module():
    from quantforge.nig import _nig_psi
    psi = levy_psi("nig", (15.0, -5.0, 0.5))
    for u in (0.5, 1.0, 2.0):
        assert psi(u) == pytest.approx(_nig_psi(u, 15.0, -5.0, 0.5))


def test_unknown_model_raises():
    with pytest.raises(ValueError):
        LevySurface("heston", (1, 2, 3), S, R, Q)
