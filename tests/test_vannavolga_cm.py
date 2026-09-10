"""Castagna-Mercurio first/second-order vanna-volga vol approximation."""

import pytest

from quantforge import VannaVolgaSmile


S, T, RD, RF = 1.30, 1.0, 0.02, 0.01
ATM, RR, BF = 0.10, -0.015, 0.004


def _vv():
    return VannaVolgaSmile(S, T, RD, RF, ATM, RR, BF)


@pytest.mark.parametrize("order", [1, 2])
def test_exact_at_pillars(order):
    vv = _vv()
    for K, sig in vv.pillars():
        assert vv.vol_cm(K, order=order) == pytest.approx(sig, abs=1e-6)


def test_second_order_matches_price_corrected_vol():
    vv = _vv()
    for K in (1.25, 1.30, 1.35):
        assert vv.vol_cm(K, order=2) == pytest.approx(
            vv.vol_price_corrected(K, RD, RF), abs=1e-3)


def test_orders_differ_off_pillar():
    vv = _vv()
    _, ka, kc = vv._ks
    mid = 0.5 * (ka + kc)
    assert abs(vv.vol_cm(mid, 1) - vv.vol_cm(mid, 2)) > 1e-6


def test_negative_rr_downward_skew():
    vv = _vv()
    kp, _, kc = vv._ks
    assert vv.vol_cm(kp, 2) > vv.vol_cm(kc, 2)


def test_first_order_is_weighted_average():
    # First order = sigma_atm + sum x_i (sigma_i - sigma_atm); at the ATM pillar
    # the weights collapse so it returns the ATM vol.
    vv = _vv()
    _, ka, _ = vv._ks
    assert vv.vol_cm(ka, 1) == pytest.approx(ATM, abs=1e-9)
