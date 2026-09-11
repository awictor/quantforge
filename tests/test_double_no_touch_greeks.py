"""Greeks of a double-no-touch (double_no_touch_greeks)."""

import pytest

from quantforge import double_no_touch_greeks, double_no_touch


S, T, R = 100.0, 0.5, 0.03
L, U, SIG = 88.0, 116.0, 0.2


def test_delta_matches_independent_bump():
    g = double_no_touch_greeks(S, L, U, T, R, SIG)
    h = 0.05
    up = double_no_touch(S + h, L, U, T, R, SIG)
    dn = double_no_touch(S - h, L, U, T, R, SIG)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-4)


def test_vega_negative():
    # A DNT is short volatility: more vol raises the knock probability.
    g = double_no_touch_greeks(S, L, U, T, R, SIG)
    assert g["vega"] < 0.0


def test_gamma_negative_inside_band():
    g = double_no_touch_greeks(S, L, U, T, R, SIG)
    assert g["gamma"] < 0.0


def test_theta_positive():
    # Less remaining time -> less chance to knock -> value rises, so -dV/dt > 0.
    g = double_no_touch_greeks(S, L, U, T, R, SIG)
    assert g["theta"] > 0.0


def test_barrier_sensitivity_signs():
    # Widening the band raises the value: raising L (narrowing) lowers it
    # (dV/dL < 0); raising U (widening) raises it (dV/dU > 0).
    g = double_no_touch_greeks(S, L, U, T, R, SIG)
    assert g["dV_dL"] < 0.0
    assert g["dV_dU"] > 0.0


def test_barrier_sensitivities_match_bumps():
    g = double_no_touch_greeks(S, L, U, T, R, SIG)
    hL = 0.05
    dL = (double_no_touch(S, L + hL, U, T, R, SIG)
          - double_no_touch(S, L - hL, U, T, R, SIG)) / (2 * hL)
    assert g["dV_dL"] == pytest.approx(dL, abs=1e-3)
    hU = 0.05
    dU = (double_no_touch(S, L, U + hU, T, R, SIG)
          - double_no_touch(S, L, U - hU, T, R, SIG)) / (2 * hU)
    assert g["dV_dU"] == pytest.approx(dU, abs=1e-3)


def test_price_field_matches_pricer():
    g = double_no_touch_greeks(S, L, U, T, R, SIG)
    assert g["price"] == pytest.approx(double_no_touch(S, L, U, T, R, SIG),
                                       abs=1e-12)


def test_requires_ordering():
    with pytest.raises(ValueError):
        double_no_touch_greeks(S, 120.0, 90.0, T, R, SIG)
