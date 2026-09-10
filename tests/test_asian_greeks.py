"""Tests for Asian option Greeks (finite differences)."""

import pytest

from quantforge import asian_greeks, geometric_asian, arithmetic_asian, OptionType


def test_geometric_delta_matches_re_difference():
    S, K = 100, 100
    g = asian_greeks(S, K, 1.0, 0.05, 0.3, OptionType.CALL, average="geometric")
    h = 0.05
    up = geometric_asian(S + h, K, 1.0, 0.05, 0.3, OptionType.CALL)
    dn = geometric_asian(S - h, K, 1.0, 0.05, 0.3, OptionType.CALL)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-3)


def test_arithmetic_delta_matches_re_difference():
    S, K = 100, 100
    g = asian_greeks(S, K, 1.0, 0.05, 0.3, OptionType.CALL, average="arithmetic")
    h = 0.05
    up = arithmetic_asian(S + h, K, 1.0, 0.05, 0.3, OptionType.CALL)
    dn = arithmetic_asian(S - h, K, 1.0, 0.05, 0.3, OptionType.CALL)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-3)


def test_gamma_and_vega_positive():
    for avg in ("geometric", "arithmetic"):
        g = asian_greeks(100, 100, 1.0, 0.05, 0.3, OptionType.CALL, average=avg)
        assert g["gamma"] > 0
        assert g["vega"] > 0


def test_arithmetic_delta_above_geometric():
    # The arithmetic-average call is worth more, and its delta is larger too.
    geo = asian_greeks(100, 100, 1.0, 0.05, 0.3, OptionType.CALL, average="geometric")
    ari = asian_greeks(100, 100, 1.0, 0.05, 0.3, OptionType.CALL, average="arithmetic")
    assert ari["delta"] > geo["delta"]


def test_call_delta_positive_put_delta_negative():
    c = asian_greeks(100, 100, 1.0, 0.05, 0.25, OptionType.CALL)
    p = asian_greeks(100, 100, 1.0, 0.05, 0.25, OptionType.PUT)
    assert c["delta"] > 0
    assert p["delta"] < 0


def test_price_field_matches_direct():
    g = asian_greeks(100, 95, 0.5, 0.04, 0.3, OptionType.CALL, average="geometric")
    assert g["price"] == pytest.approx(
        geometric_asian(100, 95, 0.5, 0.04, 0.3, OptionType.CALL), abs=1e-9)


def test_rejects_unknown_average():
    with pytest.raises(ValueError):
        asian_greeks(100, 100, 1.0, 0.05, 0.3, average="harmonic")
