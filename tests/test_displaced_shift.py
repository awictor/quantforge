"""Tests for displaced-diffusion shift calibration from a vol skew."""

import math

import pytest

from quantforge import (
    OptionType,
    displaced_diffusion_price,
    displaced_implied_shift,
    implied_volatility,
)


S, T, R = 100.0, 1.0, 0.05
STRIKES = [80, 90, 100, 110, 120]


def _quotes_from_shift(shift, sigma=0.2):
    out = []
    for K in STRIKES:
        px = displaced_diffusion_price(S, K, T, R, sigma, shift,
                                       OptionType.CALL, b=R)
        iv = implied_volatility(px, S, K, T, R, OptionType.CALL, b=R)
        out.append((K, iv))
    return out


@pytest.mark.parametrize("true_shift", [0.0, 40.0, 100.0, 200.0])
def test_recovers_known_shift(true_shift):
    quotes = _quotes_from_shift(true_shift)
    shift, sigma_atm, rmse = displaced_implied_shift(S, T, R, quotes, b=R)
    assert shift == pytest.approx(true_shift, abs=0.5)
    assert sigma_atm == pytest.approx(0.2, abs=1e-4)
    assert rmse < 1e-4


def test_positive_shift_for_downward_skew():
    # A downward skew (low-strike vol > high-strike vol) implies a positive shift.
    quotes = _quotes_from_shift(60.0)
    lo_iv = quotes[0][1]
    hi_iv = quotes[-1][1]
    assert lo_iv > hi_iv
    shift, _, _ = displaced_implied_shift(S, T, R, quotes, b=R)
    assert shift > 0


def test_flat_smile_gives_zero_shift():
    quotes = [(K, 0.2) for K in STRIKES]
    shift, sigma_atm, rmse = displaced_implied_shift(S, T, R, quotes, b=R)
    assert shift == pytest.approx(0.0, abs=0.5)
    assert sigma_atm == pytest.approx(0.2, abs=1e-4)
    assert rmse < 1e-6


def test_requires_two_quotes():
    with pytest.raises(ValueError):
        displaced_implied_shift(S, T, R, [(100, 0.2)], b=R)


def test_rejects_nonpositive_vol():
    with pytest.raises(ValueError):
        displaced_implied_shift(S, T, R, [(90, 0.2), (110, 0.0)], b=R)


def test_refit_reproduces_quotes():
    # Round-trip: calibrate, then re-price at the fitted (shift, sigma) and
    # confirm the model implied vols match the input quotes.
    quotes = _quotes_from_shift(80.0)
    shift, sigma_atm, _ = displaced_implied_shift(S, T, R, quotes, b=R)
    for K, iv_mkt in quotes:
        px = displaced_diffusion_price(S, K, T, R, sigma_atm, shift,
                                       OptionType.CALL, b=R)
        iv_mod = implied_volatility(px, S, K, T, R, OptionType.CALL, b=R)
        assert iv_mod == pytest.approx(iv_mkt, abs=1e-4)
