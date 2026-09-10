"""Tests for Bachelier caps, floors, and collars."""

import math

import pytest

from quantforge import (
    CapletPeriod, caplet_price, cap_price, floor_price, collar_price,
    caplet_floorlet_parity,
)


def _strip():
    # A simple 4-period annual strip: rising forwards, flat vol and discounting.
    periods = []
    for i in range(1, 5):
        df = math.exp(-0.03 * i)          # discount to payment date
        periods.append(CapletPeriod(forward=0.03 + 0.002 * i, expiry=float(i),
                                    accrual=1.0, discount=df, sigma_n=0.01))
    return periods


def test_caplet_floorlet_parity():
    p = _strip()[0]
    K = 0.035
    caplet = caplet_price(p, K, is_cap=True)
    floorlet = caplet_price(p, K, is_cap=False)
    assert caplet - floorlet == pytest.approx(caplet_floorlet_parity(p, K), abs=1e-12)


def test_cap_minus_floor_equals_swap_pv():
    # Cap(K) - Floor(K) = sum of discounted accrual*(F - K) over the strip.
    periods = _strip()
    K = 0.035
    cap = cap_price(periods, K)
    floor = floor_price(periods, K)
    swap_pv = sum(p.discount * p.accrual * (p.forward - K) for p in periods)
    assert cap - floor == pytest.approx(swap_pv, abs=1e-10)


def test_cap_positive_and_increasing_in_vol():
    periods = _strip()
    lo = cap_price(periods, 0.04)
    hi_vol = [CapletPeriod(p.forward, p.expiry, p.accrual, p.discount, p.sigma_n * 3)
              for p in periods]
    assert lo > 0
    assert cap_price(hi_vol, 0.04) > lo


def test_cap_decreasing_in_strike():
    periods = _strip()
    assert cap_price(periods, 0.03) > cap_price(periods, 0.05)


def test_floor_increasing_in_strike():
    periods = _strip()
    assert floor_price(periods, 0.05) > floor_price(periods, 0.03)


def test_handles_negative_forward_rates():
    # Bachelier permits negative rates; a floor struck at 0 still has value.
    periods = [CapletPeriod(forward=-0.005, expiry=1.0, accrual=1.0,
                            discount=math.exp(-0.0), sigma_n=0.01)]
    fl = floor_price(periods, 0.0)
    assert fl > 0


def test_collar_is_cap_minus_floor():
    periods = _strip()
    c = collar_price(periods, cap_strike=0.045, floor_strike=0.03)
    assert c == pytest.approx(cap_price(periods, 0.045) - floor_price(periods, 0.03),
                              abs=1e-12)


def test_zero_cost_collar_exists_between_strikes():
    # Collar value moves monotonically from positive to negative as we raise the
    # cap strike, so a zero-cost strike sits in between.
    periods = _strip()
    low = collar_price(periods, cap_strike=0.03, floor_strike=0.03)   # cap rich
    high = collar_price(periods, cap_strike=0.06, floor_strike=0.03)  # cap cheap
    assert low > 0 > high
