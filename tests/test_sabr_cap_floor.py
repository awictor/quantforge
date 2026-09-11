"""SABR-smile cap / floor pricing (rates.sabr_cap_price, sabr_floor_price)."""

import pytest

from quantforge import (
    sabr_cap_price, sabr_floor_price, caplet_price, cap_price, CapletPeriod,
)
from quantforge.sabr import sabr_normal_vol


def _periods():
    return [CapletPeriod(forward=0.03 + 0.002 * i, expiry=i + 1, accrual=1.0,
                         discount=1.0 / (1.05 ** (i + 1)), sigma_n=0.01)
            for i in range(5)]


K = 0.04
ALPHA, BETA, RHO, NU = 0.01, 0.0, -0.2, 0.3


def test_cap_equals_sum_of_sabr_caplets():
    periods = _periods()
    manual = 0.0
    for p in periods:
        v = sabr_normal_vol(p.forward, K, p.expiry, ALPHA, BETA, RHO, NU)
        pp = CapletPeriod(forward=p.forward, expiry=p.expiry, accrual=p.accrual,
                          discount=p.discount, sigma_n=v)
        manual += caplet_price(pp, K, is_cap=True)
    assert sabr_cap_price(periods, K, ALPHA, BETA, RHO, NU) == pytest.approx(
        manual, abs=1e-12)


def test_cap_floor_parity():
    periods = _periods()
    cap = sabr_cap_price(periods, K, ALPHA, BETA, RHO, NU)
    flr = sabr_floor_price(periods, K, ALPHA, BETA, RHO, NU)
    parity = sum(p.discount * p.accrual * (p.forward - K) for p in periods)
    assert cap - flr == pytest.approx(parity, abs=1e-10)


def test_flat_sabr_close_to_normal_cap():
    # beta = nu = 0: SABR normal vol is ~alpha, so the cap is close to a flat
    # normal cap at sigma_n = alpha (exact only at the money per caplet).
    periods = _periods()
    flat = [CapletPeriod(forward=p.forward, expiry=p.expiry, accrual=p.accrual,
                         discount=p.discount, sigma_n=ALPHA) for p in periods]
    sabr = sabr_cap_price(periods, K, ALPHA, 0.0, 0.0, 0.0)
    normal = cap_price(flat, K)
    assert sabr == pytest.approx(normal, rel=0.02)


def test_prices_positive():
    periods = _periods()
    assert sabr_cap_price(periods, K, ALPHA, BETA, RHO, NU) > 0.0
    assert sabr_floor_price(periods, K, ALPHA, BETA, RHO, NU) > 0.0


def test_zero_vol_of_vol_does_not_crash():
    # nu = 0 previously divided by zero in the normal-vol expansion.
    periods = _periods()
    v = sabr_cap_price(periods, K, ALPHA, 0.5, 0.0, 0.0)
    assert v > 0.0
