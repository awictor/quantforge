"""SABR-smile European swaption pricing (rates.sabr_swaption_price)."""

import pytest

from quantforge import (
    sabr_swaption_price, swaption_implied_black_vol,
    swaption_implied_normal_vol, sabr_vol, sabr_normal_vol, CapletPeriod,
)


def _periods():
    return [CapletPeriod(forward=0.04, expiry=i + 1, accrual=1.0,
                         discount=1.0 / (1.05 ** (i + 1)), sigma_n=0.01)
            for i in range(5)]


SR, T = 0.04, 2.0
ALPHA, BETA, RHO, NU = 0.2, 0.5, -0.3, 0.4


@pytest.mark.parametrize("K", [0.03, 0.04, 0.05])
def test_black_price_recovers_sabr_vol(K):
    p = _periods()
    price = sabr_swaption_price(SR, K, T, p, ALPHA, BETA, RHO, NU, True, "black")
    iv = swaption_implied_black_vol(price, SR, K, T, p, True)
    assert iv == pytest.approx(sabr_vol(SR, K, T, ALPHA, BETA, RHO, NU), abs=1e-6)


@pytest.mark.parametrize("K", [0.03, 0.04, 0.05])
def test_normal_price_recovers_sabr_normal_vol(K):
    p = _periods()
    price = sabr_swaption_price(SR, K, T, p, ALPHA, BETA, RHO, NU, True, "normal")
    iv = swaption_implied_normal_vol(price, SR, K, T, p, True)
    assert iv == pytest.approx(sabr_normal_vol(SR, K, T, ALPHA, BETA, RHO, NU),
                               abs=1e-6)


def test_prices_positive():
    p = _periods()
    for model in ("black", "normal"):
        assert sabr_swaption_price(SR, 0.04, T, p, ALPHA, BETA, RHO, NU, True,
                                   model) > 0.0


def test_payer_above_receiver_when_itm():
    # Payer struck below the forward is in the money vs the receiver.
    p = _periods()
    pay = sabr_swaption_price(SR, 0.03, T, p, ALPHA, BETA, RHO, NU, True, "black")
    rec = sabr_swaption_price(SR, 0.03, T, p, ALPHA, BETA, RHO, NU, False, "black")
    assert pay > rec


def test_bad_model_raises():
    with pytest.raises(ValueError):
        sabr_swaption_price(SR, 0.04, T, _periods(), ALPHA, BETA, RHO, NU, True, "x")
