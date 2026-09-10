"""Greeks of a G2++ zero-coupon-bond option (g2pp_bond_option_greeks)."""

import pytest

from quantforge import g2pp_bond_option_greeks, g2pp_bond_option


P0S, P0T = 0.97, 0.94
A, B, SIG, ETA, RHO = 0.1, 0.05, 0.01, 0.008, -0.7
EXP, MAT, K = 1.0, 3.0, 0.96


def test_delta_T_matches_finite_difference():
    g = g2pp_bond_option_greeks(P0S, P0T, A, B, SIG, ETA, RHO, EXP, MAT, K, True)
    h = 1e-6
    fd = (g2pp_bond_option(P0S, P0T + h, A, B, SIG, ETA, RHO, EXP, MAT, K)
          - g2pp_bond_option(P0S, P0T - h, A, B, SIG, ETA, RHO, EXP, MAT, K)) / (2 * h)
    assert g["delta_T"] == pytest.approx(fd, abs=1e-6)


def test_delta_S_matches_finite_difference():
    g = g2pp_bond_option_greeks(P0S, P0T, A, B, SIG, ETA, RHO, EXP, MAT, K, True)
    h = 1e-6
    fd = (g2pp_bond_option(P0S + h, P0T, A, B, SIG, ETA, RHO, EXP, MAT, K)
          - g2pp_bond_option(P0S - h, P0T, A, B, SIG, ETA, RHO, EXP, MAT, K)) / (2 * h)
    assert g["delta_S"] == pytest.approx(fd, abs=1e-5)


def test_call_delta_T_in_unit_interval():
    g = g2pp_bond_option_greeks(P0S, P0T, A, B, SIG, ETA, RHO, EXP, MAT, K, True)
    assert 0.0 < g["delta_T"] < 1.0
    assert g["delta_S"] < 0.0


def test_vegas_positive():
    g = g2pp_bond_option_greeks(P0S, P0T, A, B, SIG, ETA, RHO, EXP, MAT, K, True)
    assert g["vega_sigma"] > 0.0
    assert g["vega_eta"] > 0.0


def test_price_field_matches_bond_option():
    g = g2pp_bond_option_greeks(P0S, P0T, A, B, SIG, ETA, RHO, EXP, MAT, K, True)
    assert g["price"] == pytest.approx(
        g2pp_bond_option(P0S, P0T, A, B, SIG, ETA, RHO, EXP, MAT, K, True),
        abs=1e-12)


def test_put_delta_T_negative():
    g = g2pp_bond_option_greeks(P0S, P0T, A, B, SIG, ETA, RHO, EXP, MAT, K, False)
    assert g["delta_T"] < 0.0
