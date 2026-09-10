"""Tests for the Bachelier (normal) model."""

import math

import pytest

from quantforge import (
    bachelier_price, bachelier_delta, bachelier_gamma, bachelier_vega,
    bachelier_implied_vol, OptionType,
)


def test_atm_closed_form():
    # At F == K with r=0 the call value is sigma * sqrt(t / (2 pi)).
    F = K = 100.0
    t, sigma = 1.0, 5.0
    expected = sigma * math.sqrt(t / (2.0 * math.pi))
    assert bachelier_price(F, K, t, 0.0, sigma, OptionType.CALL) == pytest.approx(expected, abs=1e-10)


def test_put_call_parity():
    # C - P = e^{-rt} (F - K).
    F, K, t, r, sigma = 105.0, 100.0, 0.5, 0.03, 8.0
    c = bachelier_price(F, K, t, r, sigma, OptionType.CALL)
    p = bachelier_price(F, K, t, r, sigma, OptionType.PUT)
    assert c - p == pytest.approx(math.exp(-r * t) * (F - K), abs=1e-10)


def test_handles_negative_forward_and_strike():
    # Bachelier is defined for negative rates/spreads (lognormal is not).
    c = bachelier_price(-0.5, -1.0, 1.0, 0.0, 2.0, OptionType.CALL)
    p = bachelier_price(-0.5, -1.0, 1.0, 0.0, 2.0, OptionType.PUT)
    assert c > 0 and p > 0
    assert c - p == pytest.approx(-0.5 - (-1.0), abs=1e-10)


def test_delta_matches_fd():
    F, K, t, r, sigma = 100.0, 105.0, 0.75, 0.02, 6.0
    h = 1e-4
    for ot in (OptionType.CALL, OptionType.PUT):
        fd = (bachelier_price(F + h, K, t, r, sigma, ot)
              - bachelier_price(F - h, K, t, r, sigma, ot)) / (2 * h)
        assert bachelier_delta(F, K, t, r, sigma, ot) == pytest.approx(fd, abs=1e-6)


def test_gamma_matches_fd():
    F, K, t, r, sigma = 100.0, 105.0, 0.75, 0.02, 6.0
    h = 1e-2
    up = bachelier_price(F + h, K, t, r, sigma, OptionType.CALL)
    mid = bachelier_price(F, K, t, r, sigma, OptionType.CALL)
    dn = bachelier_price(F - h, K, t, r, sigma, OptionType.CALL)
    fd = (up - 2 * mid + dn) / (h * h)
    assert bachelier_gamma(F, K, t, r, sigma) == pytest.approx(fd, abs=1e-5)


def test_vega_matches_fd():
    F, K, t, r, sigma = 100.0, 105.0, 0.75, 0.02, 6.0
    h = 1e-4
    fd = (bachelier_price(F, K, t, r, sigma + h, OptionType.CALL)
          - bachelier_price(F, K, t, r, sigma - h, OptionType.CALL)) / (2 * h)
    assert bachelier_vega(F, K, t, r, sigma) == pytest.approx(fd, abs=1e-6)


@pytest.mark.parametrize("sigma", [1.0, 5.0, 20.0])
@pytest.mark.parametrize("K", [80, 100, 120])
@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_implied_vol_roundtrip(sigma, K, ot):
    F, t, r = 100.0, 1.0, 0.02
    price = bachelier_price(F, K, t, r, sigma, ot)
    iv = bachelier_implied_vol(price, F, K, t, r, ot)
    # Deep OTM/ITM at low normal vol carries essentially zero vega, so vol is
    # numerically unidentifiable there; assert the recovered *price* matches.
    if bachelier_vega(F, K, t, r, sigma) > 1e-8:
        assert iv == pytest.approx(sigma, rel=1e-4)
    else:
        assert bachelier_price(F, K, t, r, iv, ot) == pytest.approx(price, abs=1e-8)


def test_implied_vol_rejects_below_intrinsic():
    with pytest.raises(ValueError):
        bachelier_implied_vol(0.01, 120, 100, 1.0, 0.0, OptionType.CALL)  # < 20 intrinsic


def test_zero_time_is_intrinsic():
    assert bachelier_price(110, 100, 0.0, 0.05, 5.0, OptionType.CALL) == pytest.approx(10.0)
    assert bachelier_price(90, 100, 0.0, 0.05, 5.0, OptionType.PUT) == pytest.approx(10.0)
