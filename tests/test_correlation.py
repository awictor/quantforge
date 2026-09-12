"""Tests for index implied correlation."""

import math

import pytest

from quantforge import (
    implied_correlation, index_vol_from_correlation, dispersion_basket_vol,
    dispersion_trade_pnl,
)

W = [0.4, 0.35, 0.25]
V = [0.25, 0.30, 0.20]


def test_rho_one_is_weighted_average_vol():
    # Perfectly correlated members: index vol = sum(w_i sigma_i).
    idx = index_vol_from_correlation(W, V, 1.0)
    assert idx == pytest.approx(sum(w * s for w, s in zip(W, V)), abs=1e-9)


def test_rho_zero_is_dispersion_basket_vol():
    idx = index_vol_from_correlation(W, V, 0.0)
    assert idx == pytest.approx(dispersion_basket_vol(W, V), abs=1e-9)


@pytest.mark.parametrize("rho", [-0.3, 0.0, 0.25, 0.6, 1.0])
def test_round_trip(rho):
    idx = index_vol_from_correlation(W, V, rho)
    assert implied_correlation(W, V, idx) == pytest.approx(rho, abs=1e-9)


def test_higher_index_vol_higher_correlation():
    lo = implied_correlation(W, V, 0.17)
    hi = implied_correlation(W, V, 0.24)
    assert hi > lo


def test_implied_in_unit_band_for_reasonable_quote():
    # An index vol between the dispersion (rho=0) and weighted-avg (rho=1) vols
    # implies a correlation in [0, 1].
    disp = dispersion_basket_vol(W, V)
    wavg = sum(w * s for w, s in zip(W, V))
    mid = 0.5 * (disp + wavg)
    rho = implied_correlation(W, V, mid)
    assert 0.0 <= rho <= 1.0


def test_two_name_basket():
    w, v = [0.5, 0.5], [0.2, 0.2]
    # Equal names, equal vols: index vol^2 = 0.5*(1+rho)*sigma^2.
    idx = index_vol_from_correlation(w, v, 0.5)
    expected = math.sqrt(0.5 * (1 + 0.5) * 0.2 * 0.2)
    assert idx == pytest.approx(expected, abs=1e-9)


def test_rejects_bad_inputs():
    with pytest.raises(ValueError):
        implied_correlation([1.0], [0.2], 0.2)          # need >= 2 names
    with pytest.raises(ValueError):
        index_vol_from_correlation(W, V, 1.5)           # rho out of range
    with pytest.raises(ValueError):
        implied_correlation(W, V, -0.1)                 # negative index vol


DW = [0.4, 0.35, 0.25]
DMV = [0.3, 0.25, 0.35]


def test_dispersion_profits_when_realized_correlation_low():
    k_idx = index_vol_from_correlation(DW, DMV, 0.5)
    realized_idx = index_vol_from_correlation(DW, DMV, 0.3)  # calmer index
    assert dispersion_trade_pnl(DW, DMV, realized_idx, DMV, k_idx) > 0


def test_dispersion_loses_when_realized_correlation_high():
    k_idx = index_vol_from_correlation(DW, DMV, 0.5)
    realized_idx = index_vol_from_correlation(DW, DMV, 0.8)
    assert dispersion_trade_pnl(DW, DMV, realized_idx, DMV, k_idx) < 0


def test_dispersion_zero_at_strikes():
    k_idx = index_vol_from_correlation(DW, DMV, 0.5)
    assert dispersion_trade_pnl(DW, DMV, k_idx, DMV, k_idx) == pytest.approx(0.0, abs=1e-12)


def test_dispersion_validation():
    k_idx = index_vol_from_correlation(DW, DMV, 0.5)
    with pytest.raises(ValueError):
        dispersion_trade_pnl(DW, DMV, 0.2, [0.3, 0.25], k_idx)
