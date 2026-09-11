"""Black (lognormal) European swaption (rates.black_swaption_price)."""

import pytest

from quantforge import (
    black_swaption_price, black_swaption_greeks, swaption_price,
    swaption_parity, CapletPeriod, annuity,
)


def _periods():
    return [CapletPeriod(forward=0.04, expiry=i + 1, accrual=1.0,
                         discount=1.0 / (1.05 ** (i + 1)), sigma_n=0.01)
            for i in range(5)]


SR, K, T, SIGB = 0.04, 0.04, 2.0, 0.25


def test_payer_receiver_parity():
    p = _periods()
    pay = black_swaption_price(SR, K, T, SIGB, p, True)
    rec = black_swaption_price(SR, K, T, SIGB, p, False)
    assert pay - rec == pytest.approx(swaption_parity(SR, K, p), abs=1e-10)


def test_atm_close_to_bachelier():
    # ATM Black with sigma_b ~ Bachelier with sigma_n = sigma_b * swap_rate.
    p = _periods()
    black = black_swaption_price(SR, K, T, SIGB, p, True)
    bach = swaption_price(SR, K, T, SIGB * SR, p, True)
    assert black == pytest.approx(bach, rel=0.02)


def test_rate_delta_matches_fd():
    p = _periods()
    g = black_swaption_greeks(SR, K, T, SIGB, p, True)
    h = 1e-6
    fd = (black_swaption_price(SR + h, K, T, SIGB, p, True)
          - black_swaption_price(SR - h, K, T, SIGB, p, True)) / (2 * h)
    assert g["rate_delta"] == pytest.approx(fd, abs=1e-4)


def test_rate_gamma_and_vega_positive():
    g = black_swaption_greeks(SR, K, T, SIGB, _periods(), True)
    assert g["rate_gamma"] > 0.0
    assert g["vega"] > 0.0
    assert g["annuity"] == pytest.approx(annuity(_periods()), abs=1e-12)


def test_payer_delta_positive_receiver_negative():
    p = _periods()
    assert black_swaption_greeks(SR, K, T, SIGB, p, True)["rate_delta"] > 0.0
    assert black_swaption_greeks(SR, K, T, SIGB, p, False)["rate_delta"] < 0.0


def test_expiry_zero_is_intrinsic():
    p = _periods()
    itm = black_swaption_price(0.05, K, 0.0, SIGB, p, True)
    assert itm == pytest.approx(annuity(p) * (0.05 - K), abs=1e-12)


def test_negative_rate_raises():
    with pytest.raises(ValueError):
        black_swaption_price(-0.01, K, T, SIGB, _periods(), True)
    with pytest.raises(ValueError):
        black_swaption_price(SR, -0.01, T, SIGB, _periods(), True)
