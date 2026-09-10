"""Smoke test for the end-to-end vol-surface example.

Ensures the full pipeline (implied vol -> per-expiry SVI fit -> surface ->
calendar check -> local vol) runs, fits each slice well, and produces a
calendar-arbitrage-free surface on the synthetic chain.
"""

import sys
import os
import math

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))

import vol_surface as ex  # noqa: E402

from quantforge import (
    implied_volatility, calibrate_svi, VolSurface, SurfaceSlice, OptionType,
)


def _build_surface():
    chain = ex.build_chain()
    slices = []
    for T, rows in chain.items():
        ks, tv = [], []
        for K, price in rows:
            iv = implied_volatility(price, ex.SPOT, K, T, ex.RATE, OptionType.CALL)
            ks.append(math.log(K / ex.SPOT))
            tv.append(iv * iv * T)
        params, rmse = calibrate_svi(ks, tv)
        slices.append((T, SurfaceSlice(t=T, params=params, rmse=rmse), rmse))
    return slices


@pytest.fixture(scope="module")
def surface_slices():
    # SVI calibration is the expensive step; build it once for the module.
    return _build_surface()


def test_every_slice_fits_well(surface_slices):
    for T, _slice, rmse in surface_slices:
        # Multi-start calibration should fit the smooth synthetic smile tightly.
        assert rmse < 1e-4, f"expiry {T} fit poorly: rmse={rmse}"


def test_surface_is_calendar_arbitrage_free(surface_slices):
    surface = VolSurface([s for _, s, _ in surface_slices])
    assert surface.is_calendar_arbitrage_free()


def test_example_main_runs(capsys):
    ex.main()
    out = capsys.readouterr().out
    assert "Interpolated implied vol" in out
    assert "Dupire local vol" in out
    assert "violations: 0" in out
