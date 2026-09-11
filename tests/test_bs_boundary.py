"""Bjerksund-Stensland (2002) exercise boundary (bjerksund_stensland_boundary)."""

import pytest

from quantforge import (
    bjerksund_stensland, bjerksund_stensland_boundary, baw_critical_spot,
)


def test_no_dividend_call_has_no_boundary():
    assert bjerksund_stensland_boundary(100.0, 1.0, 0.05, 0.2, "call") is None


def test_call_boundary_above_strike_hits_intrinsic():
    I = bjerksund_stensland_boundary(100.0, 1.0, 0.05, 0.2, "call", b=0.0)
    assert I > 100.0
    assert bjerksund_stensland(I, 100.0, 1.0, 0.05, 0.2, "call", b=0.0) == pytest.approx(
        I - 100.0, abs=1e-3)


def test_put_boundary_below_strike_hits_intrinsic():
    I = bjerksund_stensland_boundary(100.0, 1.0, 0.05, 0.2, "put")
    assert I < 100.0
    assert bjerksund_stensland(I, 100.0, 1.0, 0.05, 0.2, "put") == pytest.approx(
        100.0 - I, abs=1e-3)


def test_just_below_call_boundary_is_continuation():
    I = bjerksund_stensland_boundary(100.0, 1.0, 0.05, 0.2, "call", b=0.0)
    price = bjerksund_stensland(I * 0.98, 100.0, 1.0, 0.05, 0.2, "call", b=0.0)
    intrinsic = I * 0.98 - 100.0
    assert price > intrinsic  # continuation value strictly above intrinsic


def test_boundary_in_same_ballpark_as_baw():
    # Two flat-boundary approximations; expect them within ~10%.
    bs_put = bjerksund_stensland_boundary(100.0, 1.0, 0.05, 0.2, "put")
    baw_put = baw_critical_spot(100.0, 1.0, 0.05, 0.2, "put")
    assert bs_put == pytest.approx(baw_put, rel=0.1)


def test_prices_unchanged_by_refactor():
    # Guard the trigger-extraction refactor: known BS2002 values.
    assert bjerksund_stensland(100.0, 100.0, 1.0, 0.05, 0.2, "call", b=0.0) == pytest.approx(
        7.652273, abs=1e-4)
    assert bjerksund_stensland(100.0, 100.0, 1.0, 0.05, 0.2, "put") == pytest.approx(
        6.015934, abs=1e-4)
