"""Tests for realized-volatility estimators.

We synthesize an OHLC series from a geometric Brownian motion with a known
annual volatility (seeded, deterministic) and check that every estimator
recovers it within a tolerance set by the sample size. We also pin exact
algebraic values on tiny hand-built series.
"""

import math
import random

import pytest

from quantforge import (
    close_to_close, ewma_vol, parkinson, garman_klass, rogers_satchell,
    yang_zhang, vol_report,
)


def _synth_ohlc(n=2000, sigma=0.20, mu=0.05, ppy=252, seed=0, intraday_steps=8):
    """Build daily OHLC bars from an intraday GBM with annual vol ``sigma``."""
    rng = random.Random(seed)
    dt_day = 1.0 / ppy
    dt = dt_day / intraday_steps
    drift = (mu - 0.5 * sigma * sigma) * dt
    vol = sigma * math.sqrt(dt)

    s = 100.0
    opens, highs, lows, closes = [], [], [], []
    for _ in range(n):
        o = s
        hi = lo = s
        for _ in range(intraday_steps):
            s *= math.exp(drift + vol * rng.gauss(0.0, 1.0))
            hi = max(hi, s)
            lo = min(lo, s)
        opens.append(o)
        highs.append(hi)
        lows.append(lo)
        closes.append(s)
    return opens, highs, lows, closes


# --- Exact algebraic checks on tiny series ---
def test_close_to_close_constant_growth_is_zero():
    # Constant multiplicative step -> identical log returns -> zero variance.
    closes = [100 * (1.01 ** i) for i in range(10)]
    assert close_to_close(closes) == pytest.approx(0.0, abs=1e-12)


def test_parkinson_known_value():
    # One bar with H/L = e^0.02: ln(H/L)=0.02 -> var = 0.02^2/(4 ln2).
    highs = [100 * math.exp(0.02)]
    lows = [100.0]
    v = parkinson(highs, lows, periods_per_year=1)
    expected = math.sqrt(0.02 ** 2 / (4 * math.log(2)))
    assert v == pytest.approx(expected, rel=1e-12)


def test_ewma_constant_returns():
    # Identical returns r -> EWMA variance converges to r^2.
    closes = [100 * (1.005 ** i) for i in range(200)]
    r = math.log(1.005)
    assert ewma_vol(closes, lam=0.94, periods_per_year=1) == pytest.approx(abs(r), rel=1e-6)


# --- Recovery of a known volatility from synthetic data ---
def test_all_estimators_recover_known_vol():
    sigma = 0.25
    # Dense intraday sampling so the observed bar range approaches the true
    # continuous high/low; too few steps biases range estimators downward.
    opens, highs, lows, closes = _synth_ohlc(n=3000, sigma=sigma, seed=42,
                                             intraday_steps=64)
    ests = {
        "c2c": close_to_close(closes),
        "parkinson": parkinson(highs, lows),
        "garman_klass": garman_klass(opens, highs, lows, closes),
        "rogers_satchell": rogers_satchell(opens, highs, lows, closes),
        "yang_zhang": yang_zhang(opens, highs, lows, closes),
    }
    for name, v in ests.items():
        assert v == pytest.approx(sigma, rel=0.15), f"{name}={v}"


def test_range_estimators_more_stable_than_c2c():
    # Across independent samples, Parkinson's spread should be tighter than
    # close-to-close (it is a lower-variance estimator).
    sigma = 0.3
    c2c_vals, park_vals = [], []
    for seed in range(20):
        o, h, l, c = _synth_ohlc(n=200, sigma=sigma, seed=seed)
        c2c_vals.append(close_to_close(c))
        park_vals.append(parkinson(h, l))

    def spread(xs):
        m = sum(xs) / len(xs)
        return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))

    assert spread(park_vals) < spread(c2c_vals)


def test_vol_report_bundles_all():
    o, h, l, c = _synth_ohlc(n=500, sigma=0.2, seed=7)
    rep = vol_report(o, h, l, c)
    assert rep.close_to_close > 0
    assert rep.parkinson > 0
    assert rep.yang_zhang > 0
    # Sanity: all estimators land in the same ballpark.
    vals = [rep.close_to_close, rep.parkinson, rep.garman_klass,
            rep.rogers_satchell, rep.yang_zhang]
    assert max(vals) / min(vals) < 2.0


# --- Input validation ---
def test_rejects_bad_bars():
    with pytest.raises(ValueError):
        garman_klass([100], [90], [80], [95])  # high 90 < close 95


def test_rejects_negative_prices():
    with pytest.raises(ValueError):
        close_to_close([100, -5, 110])
