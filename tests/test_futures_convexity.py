"""Futures/forward convexity adjustment."""

import pytest

from quantforge import (
    ho_lee_convexity_adjustment, hull_white_convexity_adjustment,
    forward_from_futures, futures_from_forward,
    forward_curve_from_futures_strip, stub_discount_factors_from_forwards,
)


STRIP = [(0.25, 0.5, 0.05), (0.5, 0.75, 0.052), (0.75, 1.0, 0.054),
         (1.0, 1.25, 0.056)]


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


def test_strip_forwards_below_futures():
    fwds = forward_curve_from_futures_strip(STRIP, 0.012)
    assert all(fwds[i][2] < STRIP[i][2] for i in range(len(STRIP)))


def test_strip_adjustment_grows_down_curve():
    fwds = forward_curve_from_futures_strip(STRIP, 0.012)
    adj = [STRIP[i][2] - fwds[i][2] for i in range(len(STRIP))]
    assert all(adj[i] < adj[i + 1] for i in range(len(adj) - 1))


def test_strip_matches_pointwise():
    fwds = forward_curve_from_futures_strip(STRIP, 0.012)
    for i, (t1, t2, fut) in enumerate(STRIP):
        assert fwds[i][2] == pytest.approx(forward_from_futures(fut, 0.012, t1, t2))


def test_discount_factors_decreasing():
    fwds = forward_curve_from_futures_strip(STRIP, 0.012)
    df = stub_discount_factors_from_forwards(fwds)
    assert all(df[i][1] > df[i + 1][1] for i in range(len(df) - 1))
    assert all(d[1] < 1 for d in df)
    assert df[0][1] == pytest.approx(1 / (1 + fwds[0][2] * 0.25))


def test_strip_validation():
    with pytest.raises(ValueError):
        stub_discount_factors_from_forwards([(0.5, 0.5, 0.05)])


def test_validation():
    with pytest.raises(ValueError):
        ho_lee_convexity_adjustment(SIG, 2, 1)
    with pytest.raises(ValueError):
        hull_white_convexity_adjustment(SIG, -0.05, T1, T2)
