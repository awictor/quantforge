"""Tests for the term-structure vol surface and calendar-arbitrage checks."""

import math

import pytest

from quantforge import VolSurface, SurfaceSlice, SVIParams


def _flat_slice(t, variance_level):
    # A flat smile: w(k) = variance_level for all k (a=level, b=0).
    return SurfaceSlice(t=t, params=SVIParams(a=variance_level, b=0.0, rho=0.0,
                                              m=0.0, s=0.1), rmse=0.0)


def test_interpolates_total_variance_linearly_in_t():
    # Flat smiles with total variance 0.04 at t=1 and 0.09 at t=2.
    s = VolSurface([_flat_slice(1.0, 0.04), _flat_slice(2.0, 0.09)])
    # Midway in t: linear interpolation of total variance -> 0.065.
    assert s.total_variance(0.0, 1.5) == pytest.approx(0.065, abs=1e-12)


def test_implied_vol_from_total_variance():
    s = VolSurface([_flat_slice(1.0, 0.04), _flat_slice(2.0, 0.16)])
    # At t=1, w=0.04 -> vol=0.2; at t=2, w=0.16 -> vol=sqrt(0.08).
    assert s.implied_vol(0.0, 1.0) == pytest.approx(0.2, abs=1e-9)
    assert s.implied_vol(0.0, 2.0) == pytest.approx(math.sqrt(0.16 / 2.0), abs=1e-9)


def test_calendar_arbitrage_free_surface():
    # Total variance increasing in t at every k: no arbitrage.
    s = VolSurface([_flat_slice(0.5, 0.02), _flat_slice(1.0, 0.04),
                    _flat_slice(2.0, 0.09)])
    assert s.is_calendar_arbitrage_free()
    assert s.calendar_arbitrage() == []


def test_calendar_arbitrage_detected():
    # Total variance FALLS from t=1 (0.09) to t=2 (0.04): calendar arbitrage.
    s = VolSurface([_flat_slice(1.0, 0.09), _flat_slice(2.0, 0.04)])
    violations = s.calendar_arbitrage()
    assert len(violations) > 0
    v = violations[0]
    assert v.t_short == 1.0 and v.t_long == 2.0
    assert v.w_long < v.w_short
    assert not s.is_calendar_arbitrage_free()


def test_fit_from_quotes_recovers_smiles():
    # Build synthetic per-expiry SVI quotes and fit a surface.
    truth = {
        0.5: SVIParams(a=0.02, b=0.1, rho=-0.3, m=0.0, s=0.1),
        1.0: SVIParams(a=0.04, b=0.15, rho=-0.35, m=0.0, s=0.12),
    }
    ks = [-0.3, -0.15, 0.0, 0.15, 0.3]
    quotes = [(t, ks, [p.total_variance(k) for k in ks]) for t, p in truth.items()]
    surf = VolSurface.fit(quotes)
    assert len(surf.slices) == 2
    for sl in surf.slices:
        assert sl.rmse < 1e-4
    # Fitted surface reproduces the truth at the quote points.
    for t, p in truth.items():
        for k in ks:
            assert surf.total_variance(k, t) == pytest.approx(p.total_variance(k), abs=1e-3)


def test_extrapolation_beyond_last_expiry_increases_variance():
    s = VolSurface([_flat_slice(1.0, 0.04), _flat_slice(2.0, 0.09)])
    # Total variance must keep rising past the last expiry (forward variance>0).
    assert s.total_variance(0.0, 3.0) > s.total_variance(0.0, 2.0)


def test_short_end_scales_to_zero():
    s = VolSurface([_flat_slice(1.0, 0.04)])
    # Before the first expiry, total variance scales toward 0 as t->0.
    assert s.total_variance(0.0, 0.5) == pytest.approx(0.02, abs=1e-12)
