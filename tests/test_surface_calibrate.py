"""Tests for calendar-arbitrage-free surface fitting."""

import pytest

from quantforge import VolSurface, SVIParams


KS = [-0.3, -0.15, 0.0, 0.15, 0.3]


def _quotes_from(truths):
    return [(t, KS, [p.total_variance(k) for k in KS]) for t, p in truths]


def test_repairs_calendar_arbitrage():
    # 1y level above 2y level -> total variance falls with maturity (arb).
    truths = [
        (1.0, SVIParams(a=0.10, b=0.1, rho=-0.3, m=0.0, s=0.2)),
        (2.0, SVIParams(a=0.04, b=0.1, rho=-0.3, m=0.0, s=0.2)),
    ]
    q = _quotes_from(truths)
    raw = VolSurface.fit(q)
    fixed = VolSurface.fit_arbitrage_free(q)
    assert not raw.is_calendar_arbitrage_free()
    assert fixed.is_calendar_arbitrage_free()


def test_consistent_surface_left_arbitrage_free():
    truths = [
        (1.0, SVIParams(a=0.04, b=0.1, rho=-0.3, m=0.0, s=0.2)),
        (2.0, SVIParams(a=0.09, b=0.1, rho=-0.3, m=0.0, s=0.2)),
    ]
    fixed = VolSurface.fit_arbitrage_free(_quotes_from(truths))
    assert fixed.is_calendar_arbitrage_free()


def test_repair_only_lifts_variance():
    # The repaired long slice's total variance is >= the raw fit everywhere.
    truths = [
        (1.0, SVIParams(a=0.10, b=0.1, rho=-0.3, m=0.0, s=0.2)),
        (2.0, SVIParams(a=0.04, b=0.1, rho=-0.3, m=0.0, s=0.2)),
    ]
    q = _quotes_from(truths)
    raw = VolSurface.fit(q)
    fixed = VolSurface.fit_arbitrage_free(q)
    for k in KS:
        assert fixed.slices[1].params.total_variance(k) >= \
            raw.slices[1].params.total_variance(k) - 1e-9


def test_three_expiry_surface_repaired():
    truths = [
        (0.5, SVIParams(a=0.08, b=0.1, rho=-0.3, m=0.0, s=0.2)),
        (1.0, SVIParams(a=0.05, b=0.1, rho=-0.3, m=0.0, s=0.2)),   # dips
        (2.0, SVIParams(a=0.03, b=0.1, rho=-0.3, m=0.0, s=0.2)),   # dips again
    ]
    fixed = VolSurface.fit_arbitrage_free(_quotes_from(truths))
    assert fixed.is_calendar_arbitrage_free()
    assert len(fixed.slices) == 3
