"""Andreasen-Huge single-step arbitrage-free local-vol smile."""

import math

import pytest

from quantforge import (
    andreasen_huge_prices,
    andreasen_huge_smile,
    andreasen_huge_calibrate,
)


F, T = 100.0, 1.0
STRIKES = [50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]


def test_prices_are_arbitrage_free():
    lv = [0.2] * len(STRIKES)
    pr = andreasen_huge_prices(F, STRIKES, T, lv)
    # Monotone decreasing in strike (no call-spread arbitrage).
    assert all(pr[i] >= pr[i + 1] - 1e-12 for i in range(len(pr) - 1))
    # Convex in strike (no butterfly arbitrage).
    assert all(pr[i - 1] - 2 * pr[i] + pr[i + 1] >= -1e-9
               for i in range(1, len(pr) - 1))


def test_prices_arbitrage_free_for_a_skewed_local_vol():
    # Even a steep, non-flat local vol yields a valid surface (the point of AH).
    lv = [0.35 - 0.001 * (K - 50) for K in STRIKES]
    pr = andreasen_huge_prices(F, STRIKES, T, lv)
    assert all(pr[i] >= pr[i + 1] - 1e-12 for i in range(len(pr) - 1))
    assert all(pr[i - 1] - 2 * pr[i] + pr[i + 1] >= -1e-9
               for i in range(1, len(pr) - 1))


def test_calibration_reproduces_market_skew():
    mkt = [0.28 - 0.08 * math.log(K / F) for K in STRIKES]
    lv, rmse = andreasen_huge_calibrate(F, STRIKES, T, mkt, max_iter=80)
    assert rmse < 5e-3
    sm = {round(k, 6): iv for k, iv in andreasen_huge_smile(F, STRIKES, T, lv)}
    for K in (80, 100, 120):
        k = round(math.log(K / F), 6)
        assert sm[k] == pytest.approx(0.28 - 0.08 * math.log(K / F), abs=3e-3)


def test_calibration_downward_skew_shape():
    mkt = [0.28 - 0.08 * math.log(K / F) for K in STRIKES]
    lv, _ = andreasen_huge_calibrate(F, STRIKES, T, mkt)
    sm = {round(k, 6): iv for k, iv in andreasen_huge_smile(F, STRIKES, T, lv)}
    # Compare interior strikes (the Dirichlet edge strikes invert to 0).
    lo = sm[round(math.log(80 / F), 6)]
    hi = sm[round(math.log(120 / F), 6)]
    assert lo > hi  # low strikes richer


def test_smile_sorted_and_positive():
    lv = [0.2] * len(STRIKES)
    sm = andreasen_huge_smile(F, STRIKES, T, lv)
    ks = [k for k, _ in sm]
    assert ks == sorted(ks)
    assert all(iv > 0 for _, iv in sm)


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        andreasen_huge_prices(F, [90, 110], T, [0.2, 0.2])   # < 3 strikes
    with pytest.raises(ValueError):
        andreasen_huge_prices(F, STRIKES, T, [0.2] * 3)      # mismatched length
