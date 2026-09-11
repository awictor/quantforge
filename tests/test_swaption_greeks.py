"""Greeks of a normal-model European swaption (swaption_greeks)."""

import math

import pytest

from quantforge import swaption_greeks, swaption_price, CapletPeriod


def _periods():
    return [CapletPeriod(0.03, 2.0, 1.0, math.exp(-0.03 * (2 + i)), 0.008)
            for i in range(1, 6)]


SR, K, EXP, SIG = 0.03, 0.03, 2.0, 0.008


def test_rate_delta_matches_finite_difference():
    per = _periods()
    g = swaption_greeks(SR, K, EXP, SIG, per, payer=True)
    h = 1e-6
    fd = (swaption_price(SR + h, K, EXP, SIG, per, True)
          - swaption_price(SR - h, K, EXP, SIG, per, True)) / (2 * h)
    assert g["rate_delta"] == pytest.approx(fd, abs=1e-5)


def test_payer_positive_receiver_negative():
    per = _periods()
    assert swaption_greeks(SR, K, EXP, SIG, per, payer=True)["rate_delta"] > 0.0
    assert swaption_greeks(SR, K, EXP, SIG, per, payer=False)["rate_delta"] < 0.0


def test_atm_payer_delta_is_half_annuity():
    # Bachelier ATM delta is 0.5, so the payer swaption's rate delta is
    # annuity/2 at the money.
    per = _periods()
    g = swaption_greeks(SR, SR, EXP, SIG, per, payer=True)
    assert g["rate_delta"] == pytest.approx(g["annuity"] / 2.0, abs=1e-9)


def test_vega_positive():
    per = _periods()
    assert swaption_greeks(SR, K, EXP, SIG, per, payer=True)["vega"] > 0.0


def test_price_field_matches_swaption_price():
    per = _periods()
    g = swaption_greeks(SR, K, EXP, SIG, per, payer=True)
    assert g["price"] == pytest.approx(
        swaption_price(SR, K, EXP, SIG, per, True), abs=1e-12)
