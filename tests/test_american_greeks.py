"""Tests for Bjerksund-Stensland American Greeks (finite differences)."""

import pytest

from quantforge import (
    bjerksund_stensland, bjerksund_stensland_greeks, greeks, OptionType,
)


def test_no_dividend_call_greeks_match_european():
    # b = r: American call == European call, so its Greeks equal the BSM Greeks.
    g = bjerksund_stensland_greeks(100, 100, 1.0, 0.05, 0.2, OptionType.CALL)
    e = greeks(100, 100, 1.0, 0.05, 0.2, OptionType.CALL)
    assert g["delta"] == pytest.approx(e.delta, abs=1e-4)
    assert g["gamma"] == pytest.approx(e.gamma, abs=1e-4)
    assert g["vega"] == pytest.approx(e.vega, abs=1e-2)
    assert g["theta"] == pytest.approx(e.theta, abs=1e-2)
    assert g["rho"] == pytest.approx(e.rho, abs=1e-2)


def test_price_field_matches_direct():
    g = bjerksund_stensland_greeks(100, 95, 0.5, 0.06, 0.3, OptionType.PUT, b=0.03)
    assert g["price"] == pytest.approx(
        bjerksund_stensland(100, 95, 0.5, 0.06, 0.3, OptionType.PUT, b=0.03), abs=1e-9)


def test_delta_fd_consistency():
    # The reported delta must match a coarse re-difference of the price.
    S, K, t, r, sigma, b = 100, 100, 1.0, 0.05, 0.25, 0.0
    g = bjerksund_stensland_greeks(S, K, t, r, sigma, OptionType.PUT, b=b)
    h = 0.05
    up = bjerksund_stensland(S + h, K, t, r, sigma, OptionType.PUT, b=b)
    dn = bjerksund_stensland(S - h, K, t, r, sigma, OptionType.PUT, b=b)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-3)


def test_put_delta_negative_call_delta_positive():
    gp = bjerksund_stensland_greeks(100, 100, 1.0, 0.05, 0.3, OptionType.PUT)
    gc = bjerksund_stensland_greeks(100, 100, 1.0, 0.05, 0.3, OptionType.CALL, b=-0.05)
    assert gp["delta"] < 0
    assert gc["delta"] > 0


def test_gamma_and_vega_positive():
    for ot in (OptionType.CALL, OptionType.PUT):
        g = bjerksund_stensland_greeks(100, 100, 0.75, 0.04, 0.28, ot, b=0.01)
        assert g["gamma"] > 0
        assert g["vega"] > 0


def test_all_fields_present():
    g = bjerksund_stensland_greeks(100, 100, 1.0, 0.05, 0.2, OptionType.PUT)
    assert set(g) == {"price", "delta", "gamma", "vega", "theta", "rho"}
