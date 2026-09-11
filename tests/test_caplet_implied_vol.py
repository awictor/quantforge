"""Caplet/floorlet normal implied vol (rates.caplet_implied_normal_vol)."""

import pytest

from quantforge import caplet_price, caplet_implied_normal_vol, CapletPeriod


@pytest.mark.parametrize("sn", [0.008, 0.015])
@pytest.mark.parametrize("K", [0.02, 0.03, 0.04])
@pytest.mark.parametrize("is_cap", [True, False])
def test_round_trip(sn, K, is_cap):
    p = CapletPeriod(forward=0.03, expiry=2.0, accrual=0.5, discount=0.9,
                     sigma_n=sn)
    price = caplet_price(p, K, is_cap)
    assert caplet_implied_normal_vol(price, p, K, is_cap) == pytest.approx(
        sn, abs=1e-8)


def test_negative_forward_ok():
    # Normal model handles negative rates.
    p = CapletPeriod(forward=-0.005, expiry=1.0, accrual=1.0, discount=0.99,
                     sigma_n=0.01)
    price = caplet_price(p, 0.0, is_cap=True)
    assert caplet_implied_normal_vol(price, p, 0.0, True) == pytest.approx(
        0.01, abs=1e-8)


def test_expired_raises():
    p = CapletPeriod(forward=0.03, expiry=0.0, accrual=0.5, discount=0.9,
                     sigma_n=0.01)
    with pytest.raises(ValueError):
        caplet_implied_normal_vol(0.001, p, 0.03, True)
