"""Empirical liquidity measures: Roll, Amihud, Corwin-Schultz."""

import math
import random

import pytest

from quantforge import roll_spread, amihud_illiquidity, corwin_schultz_spread


def _roll_path(spread, seed, n=5000):
    rng = random.Random(seed)
    p = 100.0
    obs = []
    for _ in range(n):
        p += 0.01 * rng.gauss(0, 1)
        q = 1 if rng.random() < 0.5 else -1
        obs.append(p + q * spread / 2)
    return obs


def test_roll_recovers_known_spread():
    spread = 0.10
    est = [roll_spread(_roll_path(spread, s)) for s in range(30)]
    assert abs(sum(est) / len(est) - spread) < 0.02


def test_roll_zero_when_no_bounce():
    # A pure random walk has no negative autocovariance -> zero implied spread.
    rng = random.Random(1)
    p = [100.0]
    for _ in range(2000):
        p.append(p[-1] + 0.01 * rng.gauss(0, 1))
    # Roll returns 0 whenever the sample autocovariance is non-negative.
    assert roll_spread(p) >= 0.0
    assert roll_spread(p) < 0.05


def test_amihud_higher_for_lower_volume():
    rets = [0.01, 0.02, 0.015]
    liquid = amihud_illiquidity(rets, [1e6, 1e6, 1e6])
    illiquid = amihud_illiquidity(rets, [1e5, 1e5, 1e5])
    assert illiquid > liquid
    assert abs(illiquid - 10 * liquid) < 1e-15    # 10x less volume -> 10x measure


def test_amihud_zero_return_zero_impact():
    assert amihud_illiquidity([0.0, 0.0], [1e6, 1e6]) == 0.0


def _cs_path(spread_pct, seed, days=500):
    rng = random.Random(seed)
    mid = 100.0
    highs, lows = [], []
    for _ in range(days):
        prices = [mid]
        for _ in range(20):
            prices.append(prices[-1] * math.exp(0.005 * rng.gauss(0, 1)))
        mid = prices[-1]
        highs.append(max(prices) * (1 + spread_pct / 2))
        lows.append(min(prices) * (1 - spread_pct / 2))
    return highs, lows


def test_corwin_schultz_positive_and_monotone():
    est = []
    for sp in (0.0, 0.005, 0.01):
        vals = [corwin_schultz_spread(*_cs_path(sp, s)) for s in range(15)]
        est.append(sum(vals) / len(vals))
    assert all(e >= 0.0 for e in est)
    assert est[0] < est[1] < est[2]        # rises with the true spread


def test_validation():
    with pytest.raises(ValueError):
        roll_spread([1.0, 2.0])
    with pytest.raises(ValueError):
        amihud_illiquidity([0.01], [0.0])
    with pytest.raises(ValueError):
        amihud_illiquidity([0.01, 0.02], [1e6])
    with pytest.raises(ValueError):
        corwin_schultz_spread([10.0], [9.0])
    with pytest.raises(ValueError):
        corwin_schultz_spread([9.0, 9.0], [10.0, 10.0])   # high < low
