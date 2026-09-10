"""Tests for the autocallable note Monte Carlo pricer."""

import math

import pytest

from quantforge import autocallable_mc


OBS = [1.0, 2.0, 3.0]


def test_never_call_no_protection_is_discounted_coupon_notional():
    # Autocall barrier unreachable, no protection: every path pays the full
    # notional plus all coupons at maturity.
    r, t, coupon = 0.03, 3.0, 0.08
    v = autocallable_mc(100, t, r, 0.25, OBS, autocall_barrier=1e9,
                        coupon=coupon, protection_barrier=None,
                        n_paths=40_000, seed=1)
    expected = math.exp(-r * t) * (1.0 + coupon * len(OBS))
    assert v.price == pytest.approx(expected, abs=1e-9)


def test_higher_coupon_raises_value():
    lo = autocallable_mc(100, 3.0, 0.03, 0.25, OBS, autocall_barrier=100,
                         coupon=0.08, protection_barrier=70, n_paths=40_000, seed=2)
    hi = autocallable_mc(100, 3.0, 0.03, 0.25, OBS, autocall_barrier=100,
                         coupon=0.15, protection_barrier=70, n_paths=40_000, seed=2)
    assert hi.price > lo.price


def test_protection_lowers_value():
    # Adding a down-and-in downside (protection barrier) reduces the value vs
    # a note that always returns the notional.
    no_prot = autocallable_mc(100, 3.0, 0.03, 0.35, OBS, autocall_barrier=120,
                              coupon=0.05, protection_barrier=None,
                              n_paths=40_000, seed=3)
    prot = autocallable_mc(100, 3.0, 0.03, 0.35, OBS, autocall_barrier=120,
                           coupon=0.05, protection_barrier=80,
                           n_paths=40_000, seed=3)
    assert prot.price < no_prot.price


def test_reproducible_with_seed():
    a = autocallable_mc(100, 3.0, 0.03, 0.25, OBS, autocall_barrier=105,
                        coupon=0.07, protection_barrier=70, n_paths=5_000, seed=9)
    b = autocallable_mc(100, 3.0, 0.03, 0.25, OBS, autocall_barrier=105,
                        coupon=0.07, protection_barrier=70, n_paths=5_000, seed=9)
    assert a.price == b.price


def test_price_positive():
    v = autocallable_mc(100, 3.0, 0.03, 0.25, OBS, autocall_barrier=100,
                        coupon=0.08, protection_barrier=70, n_paths=20_000, seed=4)
    assert v.price > 0


def test_rejects_bad_observations():
    with pytest.raises(ValueError):
        autocallable_mc(100, 3.0, 0.03, 0.25, [1.0, 1.0, 3.0],
                        autocall_barrier=100, coupon=0.08, n_paths=100)
    with pytest.raises(ValueError):
        # Last observation must equal maturity.
        autocallable_mc(100, 3.0, 0.03, 0.25, [1.0, 2.0], autocall_barrier=100,
                        coupon=0.08, n_paths=100)
