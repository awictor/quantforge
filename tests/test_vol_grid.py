"""Tests for the surface vol-grid exporters."""

import math

import pytest

from quantforge import VolSurface, SurfaceSlice, SVIParams


def _surface():
    def sl(t, w):
        return SurfaceSlice(t=t, params=SVIParams(a=w, b=0.1, rho=-0.3, m=0.0, s=0.2),
                            rmse=0.0)
    return VolSurface([sl(0.5, 0.02), sl(1.0, 0.04), sl(2.0, 0.09)])


def test_vol_grid_dimensions():
    s = _surface()
    g = s.vol_grid(ks=[-0.2, 0.0, 0.2], ts=[0.5, 1.0, 1.5])
    assert len(g["vols"]) == 3
    assert all(len(row) == 3 for row in g["vols"])
    assert g["ks"] == [-0.2, 0.0, 0.2]
    assert g["ts"] == [0.5, 1.0, 1.5]


def test_vol_grid_values_match_implied_vol():
    s = _surface()
    ks, ts = [-0.1, 0.1], [0.75, 1.25]
    g = s.vol_grid(ks, ts)
    for i, t in enumerate(ts):
        for j, k in enumerate(ks):
            assert g["vols"][i][j] == pytest.approx(s.implied_vol(k, t))


def test_strike_grid_atm_matches_surface():
    s = _surface()
    sg = s.strike_vol_grid(spot=100, strikes=[90, 100, 110], ts=[1.0])
    # Forward = spot (b=0 default), so K=100 is ATM at k=0.
    assert sg["vols"][0][1] == pytest.approx(s.implied_vol(0.0, 1.0))


def test_strike_grid_uses_carry_forward():
    s = _surface()
    b = 0.03
    sg = s.strike_vol_grid(spot=100, strikes=[100], ts=[1.0], b=b)
    fwd = 100 * math.exp(b * 1.0)
    assert sg["vols"][0][0] == pytest.approx(s.implied_vol(math.log(100 / fwd), 1.0))


def test_strike_grid_dimensions():
    s = _surface()
    sg = s.strike_vol_grid(100, [80, 90, 100, 110, 120], [0.5, 1.0])
    assert len(sg["vols"]) == 2
    assert all(len(row) == 5 for row in sg["vols"])
