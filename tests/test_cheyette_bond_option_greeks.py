"""Greeks of a Cheyette (Hull-White) bond option (cheyette_bond_option_greeks)."""

import pytest

from quantforge import cheyette_bond_option_greeks, cheyette_bond_option


P0S, P0T = 0.97, 0.90
KAPPA, SIG = 0.1, 0.01
EXP, MAT, K = 1.0, 5.0, 0.92


def test_delta_T_matches_finite_difference():
    g = cheyette_bond_option_greeks(P0S, P0T, KAPPA, SIG, EXP, MAT, K, True)
    h = 1e-6
    fd = (cheyette_bond_option(P0S, P0T + h, KAPPA, SIG, EXP, MAT, K)
          - cheyette_bond_option(P0S, P0T - h, KAPPA, SIG, EXP, MAT, K)) / (2 * h)
    assert g["delta_T"] == pytest.approx(fd, abs=1e-6)


def test_delta_S_matches_finite_difference():
    g = cheyette_bond_option_greeks(P0S, P0T, KAPPA, SIG, EXP, MAT, K, True)
    h = 1e-6
    fd = (cheyette_bond_option(P0S + h, P0T, KAPPA, SIG, EXP, MAT, K)
          - cheyette_bond_option(P0S - h, P0T, KAPPA, SIG, EXP, MAT, K)) / (2 * h)
    assert g["delta_S"] == pytest.approx(fd, abs=1e-5)


def test_call_deltas_signs():
    g = cheyette_bond_option_greeks(P0S, P0T, KAPPA, SIG, EXP, MAT, K, True)
    assert 0.0 < g["delta_T"] < 1.0
    assert g["delta_S"] < 0.0


def test_vega_positive():
    g = cheyette_bond_option_greeks(P0S, P0T, KAPPA, SIG, EXP, MAT, K, True)
    assert g["vega"] > 0.0


def test_put_delta_T_negative():
    g = cheyette_bond_option_greeks(P0S, P0T, KAPPA, SIG, EXP, MAT, K, False)
    assert g["delta_T"] < 0.0


def test_price_field_matches_bond_option():
    g = cheyette_bond_option_greeks(P0S, P0T, KAPPA, SIG, EXP, MAT, K, True)
    assert g["price"] == pytest.approx(
        cheyette_bond_option(P0S, P0T, KAPPA, SIG, EXP, MAT, K, True), abs=1e-12)
