"""Tests for lookback option Greeks (finite differences)."""

import pytest

from quantforge import (
    lookback_greeks, floating_strike_lookback, fixed_strike_lookback, OptionType,
)


def test_floating_delta_matches_re_difference():
    g = lookback_greeks(100, 1.0, 0.05, 0.3, OptionType.CALL, kind="floating")
    h = 0.05
    up = floating_strike_lookback(100 + h, 1.0, 0.05, 0.3, OptionType.CALL)
    dn = floating_strike_lookback(100 - h, 1.0, 0.05, 0.3, OptionType.CALL)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-3)


def test_fixed_delta_matches_re_difference():
    g = lookback_greeks(100, 1.0, 0.05, 0.3, OptionType.CALL, kind="fixed", K=100)
    h = 0.05
    up = fixed_strike_lookback(100 + h, 100, 1.0, 0.05, 0.3, OptionType.CALL)
    dn = fixed_strike_lookback(100 - h, 100, 1.0, 0.05, 0.3, OptionType.CALL)
    assert g["delta"] == pytest.approx((up - dn) / (2 * h), abs=1e-3)


def test_vega_positive():
    for kind, kw in (("floating", {}), ("fixed", {"K": 100})):
        g = lookback_greeks(100, 1.0, 0.05, 0.3, OptionType.CALL, kind=kind, **kw)
        assert g["vega"] > 0


def test_price_field_matches_direct():
    g = lookback_greeks(100, 0.5, 0.04, 0.3, OptionType.PUT, kind="floating")
    assert g["price"] == pytest.approx(
        floating_strike_lookback(100, 0.5, 0.04, 0.3, OptionType.PUT), abs=1e-9)


def test_floating_put_delta_negative_with_fixed_extreme():
    # Hold the running maximum fixed (a seasoned option) so the bump isolates
    # the spot sensitivity: a floating put pays S_max - S_T, so its value falls
    # as spot rises -> negative delta. (At inception s_extreme tracks spot, so
    # the total delta there is positive, which is a different quantity.)
    p = lookback_greeks(100, 1.0, 0.05, 0.3, OptionType.PUT, kind="floating",
                        s_extreme=130)
    assert p["delta"] < 0
    c = lookback_greeks(100, 1.0, 0.05, 0.3, OptionType.CALL, kind="floating",
                        s_extreme=70)
    assert c["delta"] > 0


def test_fixed_requires_strike():
    with pytest.raises(ValueError):
        lookback_greeks(100, 1.0, 0.05, 0.3, kind="fixed")


def test_rejects_unknown_kind():
    with pytest.raises(ValueError):
        lookback_greeks(100, 1.0, 0.05, 0.3, kind="asian")
