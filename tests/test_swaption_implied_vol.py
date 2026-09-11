"""Swaption implied vol: normal (Bachelier) and Black (lognormal) inversions."""

import pytest

from quantforge import (
    swaption_price, black_swaption_price,
    swaption_implied_normal_vol, swaption_implied_black_vol, CapletPeriod,
)


def _periods():
    return [CapletPeriod(forward=0.04, expiry=i + 1, accrual=1.0,
                         discount=1.0 / (1.05 ** (i + 1)), sigma_n=0.01)
            for i in range(5)]


SR, K, T = 0.04, 0.045, 2.0


@pytest.mark.parametrize("payer", [True, False])
@pytest.mark.parametrize("sig", [0.008, 0.012])
def test_normal_round_trip(payer, sig):
    p = _periods()
    price = swaption_price(SR, K, T, sig, p, payer)
    assert swaption_implied_normal_vol(price, SR, K, T, p, payer) == pytest.approx(
        sig, abs=1e-6)


@pytest.mark.parametrize("payer", [True, False])
@pytest.mark.parametrize("sig", [0.2, 0.35])
def test_black_round_trip(payer, sig):
    p = _periods()
    price = black_swaption_price(SR, K, T, sig, p, payer)
    assert swaption_implied_black_vol(price, SR, K, T, p, payer) == pytest.approx(
        sig, abs=1e-6)


def test_black_needs_positive_rate():
    with pytest.raises(ValueError):
        swaption_implied_black_vol(0.01, -0.01, K, T, _periods(), True)


def test_expired_raises():
    p = _periods()
    with pytest.raises(ValueError):
        swaption_implied_normal_vol(0.01, SR, K, 0.0, p, True)
    with pytest.raises(ValueError):
        swaption_implied_black_vol(0.01, SR, K, 0.0, p, True)
