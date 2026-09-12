"""Black <-> Bachelier volatility conversion."""

import pytest

from quantforge import black_to_normal_vol, normal_to_black_vol
from quantforge.bsm import call_price
from quantforge.bachelier import bachelier_price


def test_atm_leading_order():
    F, K, t, sb = 100, 100, 1.0, 0.2
    sn = black_to_normal_vol(F, K, t, sb)
    assert abs(sn - sb * F) < 2.0          # sigma_N ~ sigma_B * F near the money


def test_round_trip_black_normal_black():
    F, K, t, sb = 100, 100, 1.0, 0.2
    sn = black_to_normal_vol(F, K, t, sb)
    assert abs(normal_to_black_vol(F, K, t, sn) - sb) < 1e-6


def test_price_consistency():
    F, K, t, sb = 100, 100, 1.0, 0.2
    sn = black_to_normal_vol(F, K, t, sb)
    pb = call_price(F, K, t, 0.0, sb, b=0.0)
    pn = bachelier_price(F, K, t, 0.0, sn)
    assert abs(pb - pn) < 1e-6


def test_otm_and_put_round_trips():
    sn = black_to_normal_vol(100, 110, 1.0, 0.25)
    assert abs(normal_to_black_vol(100, 110, 1.0, sn) - 0.25) < 1e-6
    snp = black_to_normal_vol(100, 90, 1.0, 0.3, is_call=False)
    assert abs(normal_to_black_vol(100, 90, 1.0, snp, is_call=False) - 0.3) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        black_to_normal_vol(-1, 100, 1.0, 0.2)
    with pytest.raises(ValueError):
        normal_to_black_vol(100, 100, 1.0, -0.2)
