"""Greeks of a double-barrier knock-in call via parity (double_knock_in_call_greeks)."""

import pytest

from quantforge import (
    double_knock_in_call_greeks, double_knock_out_call_greeks,
    double_knock_in_call,
)
from quantforge.bsm import (
    delta as bs_delta, gamma as bs_gamma, vega as bs_vega,
)


S, K, T, R = 100.0, 95.0, 0.5, 0.05
L, U, SIG = 88.0, 116.0, 0.2


def test_in_out_greek_parity():
    ki = double_knock_in_call_greeks(S, K, L, U, T, R, SIG)
    ko = double_knock_out_call_greeks(S, K, L, U, T, R, SIG)
    assert ki["delta"] + ko["delta"] == pytest.approx(
        bs_delta(S, K, T, R, SIG), abs=1e-9)
    assert ki["gamma"] + ko["gamma"] == pytest.approx(
        bs_gamma(S, K, T, R, SIG), abs=1e-9)
    assert ki["vega"] + ko["vega"] == pytest.approx(
        bs_vega(S, K, T, R, SIG), abs=1e-9)


def test_barrier_sensitivities_negate_knock_out():
    ki = double_knock_in_call_greeks(S, K, L, U, T, R, SIG)
    ko = double_knock_out_call_greeks(S, K, L, U, T, R, SIG)
    assert ki["dV_dL"] == pytest.approx(-ko["dV_dL"], abs=1e-12)
    assert ki["dV_dU"] == pytest.approx(-ko["dV_dU"], abs=1e-12)


def test_barrier_sensitivity_signs():
    # Widening the corridor lowers the knock-in: raising L (narrowing) raises it
    # (dV/dL > 0), raising U (widening) lowers it (dV/dU < 0).
    ki = double_knock_in_call_greeks(S, K, L, U, T, R, SIG)
    assert ki["dV_dL"] > 0.0
    assert ki["dV_dU"] < 0.0


def test_price_field_matches_pricer():
    ki = double_knock_in_call_greeks(S, K, L, U, T, R, SIG)
    assert ki["price"] == pytest.approx(
        double_knock_in_call(S, K, L, U, T, R, SIG), abs=1e-9)


def test_requires_ordering():
    with pytest.raises(ValueError):
        double_knock_in_call_greeks(S, K, 120.0, 90.0, T, R, SIG)
