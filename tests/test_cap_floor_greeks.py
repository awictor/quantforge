"""Greeks of interest-rate caps and floors (cap_greeks / floor_greeks)."""

import math

import pytest

from quantforge import (
    cap_greeks, floor_greeks, caplet_greeks,
    cap_price, floor_price, caplet_price, CapletPeriod,
)


def _periods():
    return [CapletPeriod(forward=0.03, expiry=e, accrual=0.5,
                         discount=math.exp(-0.03 * (e + 0.5)), sigma_n=0.008)
            for e in (0.5, 1.0, 1.5, 2.0)]


K = 0.03


def test_cap_rate_delta_matches_finite_difference():
    per = _periods()
    g = cap_greeks(per, K)
    h = 1e-6
    up = [CapletPeriod(p.forward + h, p.expiry, p.accrual, p.discount, p.sigma_n)
          for p in per]
    dn = [CapletPeriod(p.forward - h, p.expiry, p.accrual, p.discount, p.sigma_n)
          for p in per]
    fd = (cap_price(up, K) - cap_price(dn, K)) / (2 * h)
    assert g["rate_delta"] == pytest.approx(fd, abs=1e-5)


def test_cap_delta_positive_floor_delta_negative():
    per = _periods()
    assert cap_greeks(per, K)["rate_delta"] > 0.0
    assert floor_greeks(per, K)["rate_delta"] < 0.0


def test_vega_positive_for_both():
    per = _periods()
    assert cap_greeks(per, K)["vega"] > 0.0
    assert floor_greeks(per, K)["vega"] > 0.0


def test_cap_price_field_matches_cap_price():
    per = _periods()
    assert cap_greeks(per, K)["price"] == pytest.approx(cap_price(per, K),
                                                        abs=1e-12)


def test_caplet_greeks_price_matches_caplet_price():
    p = _periods()[1]
    g = caplet_greeks(p, K, is_cap=True)
    assert g["price"] == pytest.approx(caplet_price(p, K, True), abs=1e-12)


def test_cap_greeks_sum_of_caplets():
    per = _periods()
    total = cap_greeks(per, K)
    summed = sum(caplet_greeks(p, K, True)["rate_delta"] for p in per)
    assert total["rate_delta"] == pytest.approx(summed, abs=1e-12)
