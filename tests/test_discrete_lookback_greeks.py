"""Greeks of a discretely-monitored fixed-strike lookback (BGK)."""

import pytest

from quantforge import (
    discrete_fixed_strike_lookback, discrete_fixed_strike_lookback_greeks,
)


S, K, R, SIG, T, N = 100.0, 100.0, 0.05, 0.2, 1.0, 50


@pytest.mark.parametrize("ot", ["call", "put"])
def test_delta_matches_fd(ot):
    g = discrete_fixed_strike_lookback_greeks(S, K, T, R, SIG, N, ot)
    h = 1e-4 * S
    fd = (discrete_fixed_strike_lookback(S + h, K, T, R, SIG, N, ot)
          - discrete_fixed_strike_lookback(S - h, K, T, R, SIG, N, ot)) / (2 * h)
    assert g["delta"] == pytest.approx(fd, abs=1e-4)


def test_price_field_matches_pricer():
    g = discrete_fixed_strike_lookback_greeks(S, K, T, R, SIG, N, "call")
    assert g["price"] == pytest.approx(
        discrete_fixed_strike_lookback(S, K, T, R, SIG, N, "call"), abs=1e-12)


def test_signs():
    c = discrete_fixed_strike_lookback_greeks(S, K, T, R, SIG, N, "call")
    p = discrete_fixed_strike_lookback_greeks(S, K, T, R, SIG, N, "put")
    assert c["delta"] > 0.0
    assert p["delta"] < 0.0
    assert c["gamma"] > 0.0
    assert c["vega"] > 0.0
    assert p["vega"] > 0.0
