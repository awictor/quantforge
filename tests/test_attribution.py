"""Tests for Greek-based P&L attribution."""

import pytest

from quantforge import attribute_pnl, delta, vega, theta, OptionType


def test_small_spot_move_residual_tiny():
    a = attribute_pnl(100, 100, 1.0, 0.05, 0.2, dS=0.5, dsigma=0.0, dt=0.0)
    assert abs(a.unexplained) < 1e-3
    assert a.explained == pytest.approx(a.delta_pnl + a.gamma_pnl + a.vega_pnl
                                        + a.theta_pnl + a.rho_pnl)


def test_delta_pnl_matches_delta_times_move():
    a = attribute_pnl(100, 100, 1.0, 0.05, 0.2, dS=1.0, dsigma=0.0, dt=0.0)
    assert a.delta_pnl == pytest.approx(delta(100, 100, 1.0, 0.05, 0.2) * 1.0)


def test_vega_pnl_matches_vega_times_dsigma():
    a = attribute_pnl(100, 100, 1.0, 0.05, 0.2, dS=0.0, dsigma=0.01, dt=0.0)
    assert a.vega_pnl == pytest.approx(vega(100, 100, 1.0, 0.05, 0.2) * 0.01)


def test_theta_pnl_matches_theta_times_dt():
    a = attribute_pnl(100, 100, 1.0, 0.05, 0.2, dS=0.0, dsigma=0.0, dt=1 / 252)
    assert a.theta_pnl == pytest.approx(theta(100, 100, 1.0, 0.05, 0.2) * (1 / 252))


def test_total_equals_explained_plus_unexplained():
    a = attribute_pnl(100, 100, 0.5, 0.05, 0.25, dS=3, dsigma=0.02, dt=1 / 252)
    assert a.total == pytest.approx(a.explained + a.unexplained, abs=1e-12)


def test_gamma_pnl_positive_for_long_option():
    # Long gamma: a move in either direction adds convexity P&L.
    up = attribute_pnl(100, 100, 0.5, 0.05, 0.25, dS=5, dsigma=0.0, dt=0.0)
    dn = attribute_pnl(100, 100, 0.5, 0.05, 0.25, dS=-5, dsigma=0.0, dt=0.0)
    assert up.gamma_pnl > 0 and dn.gamma_pnl > 0


def test_bigger_move_larger_residual():
    small = attribute_pnl(100, 100, 0.5, 0.05, 0.25, dS=2, dsigma=0.0, dt=0.0)
    big = attribute_pnl(100, 100, 0.5, 0.05, 0.25, dS=25, dsigma=0.0, dt=0.0)
    assert abs(big.unexplained) > abs(small.unexplained)


def test_quantity_scales_all_components():
    one = attribute_pnl(100, 100, 1.0, 0.05, 0.2, dS=2, dsigma=0.01, dt=0.01, qty=1)
    ten = attribute_pnl(100, 100, 1.0, 0.05, 0.2, dS=2, dsigma=0.01, dt=0.01, qty=10)
    assert ten.total == pytest.approx(10 * one.total)
    assert ten.delta_pnl == pytest.approx(10 * one.delta_pnl)
