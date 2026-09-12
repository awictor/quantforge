"""Futures/forward convexity adjustment."""

import pytest

from quantforge import (
    ho_lee_convexity_adjustment, hull_white_convexity_adjustment,
    forward_from_futures, futures_from_forward,
)


SIG, T1, T2 = 0.01, 2.0, 2.25


def test_ho_lee_formula_and_positive():
    adj = ho_lee_convexity_adjustment(SIG, T1, T2)
    assert adj > 0
    assert adj == pytest.approx(0.5 * SIG ** 2 * T1 * T2)


def test_zero_vol_no_adjustment():
    assert ho_lee_convexity_adjustment(0.0, T1, T2) == 0.0


def test_adjustment_grows_with_maturity():
    assert ho_lee_convexity_adjustment(SIG, 4, 4.25) > ho_lee_convexity_adjustment(SIG, 2, 2.25)


def test_forward_below_futures():
    assert forward_from_futures(0.05, SIG, T1, T2) < 0.05


def test_forward_futures_round_trip():
    fwd = forward_from_futures(0.05, SIG, T1, T2)
    assert futures_from_forward(fwd, SIG, T1, T2) == pytest.approx(0.05)


def test_hull_white_reduces_to_ho_lee():
    assert hull_white_convexity_adjustment(SIG, 1e-9, T1, T2) == pytest.approx(
        ho_lee_convexity_adjustment(SIG, T1, T2), abs=1e-6)


def test_mean_reversion_dampens_adjustment():
    assert hull_white_convexity_adjustment(SIG, 0.05, T1, T2) < \
        ho_lee_convexity_adjustment(SIG, T1, T2)


def test_hull_white_forward_below_futures():
    assert forward_from_futures(0.05, SIG, T1, T2, a=0.05) < 0.05


def test_validation():
    with pytest.raises(ValueError):
        ho_lee_convexity_adjustment(SIG, 2, 1)
    with pytest.raises(ValueError):
        hull_white_convexity_adjustment(SIG, -0.05, T1, T2)
