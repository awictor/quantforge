"""Greeks of an Ikeda-Kunitomo double knock-out call (double_knock_out_call_greeks)."""

import pytest

from quantforge import double_knock_out_call_greeks, double_knock_out_call


S, K, T, R = 100.0, 95.0, 0.5, 0.05
L, U, SIG = 88.0, 116.0, 0.2


def test_delta_matches_independent_bump():
    g = double_knock_out_call_greeks(S, K, L, U, T, R, SIG)
    h = 0.05
    up = double_knock_out_call(S + h, K, L, U, T, R, SIG)
    dn = double_knock_out_call(S - h, K, L, U, T, R, SIG)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-4)


def test_vega_negative():
    # A knock-out is short vol: more vol raises the barrier-hit probability.
    g = double_knock_out_call_greeks(S, K, L, U, T, R, SIG)
    assert g["vega"] < 0.0


def test_barrier_sensitivity_signs():
    # Widening the corridor raises value: dV/dL < 0 (raising L narrows it),
    # dV/dU > 0 (raising U widens it).
    g = double_knock_out_call_greeks(S, K, L, U, T, R, SIG)
    assert g["dV_dL"] < 0.0
    assert g["dV_dU"] > 0.0


def test_barrier_sensitivities_match_bumps():
    g = double_knock_out_call_greeks(S, K, L, U, T, R, SIG)
    hL = 0.05
    dL = (double_knock_out_call(S, K, L + hL, U, T, R, SIG)
          - double_knock_out_call(S, K, L - hL, U, T, R, SIG)) / (2 * hL)
    assert g["dV_dL"] == pytest.approx(dL, abs=1e-3)
    hU = 0.05
    dU = (double_knock_out_call(S, K, L, U + hU, T, R, SIG)
          - double_knock_out_call(S, K, L, U - hU, T, R, SIG)) / (2 * hU)
    assert g["dV_dU"] == pytest.approx(dU, abs=1e-3)


def test_price_field_matches_pricer():
    g = double_knock_out_call_greeks(S, K, L, U, T, R, SIG)
    assert g["price"] == pytest.approx(
        double_knock_out_call(S, K, L, U, T, R, SIG), abs=1e-12)


def test_requires_ordering():
    with pytest.raises(ValueError):
        double_knock_out_call_greeks(S, K, 120.0, 90.0, T, R, SIG)
