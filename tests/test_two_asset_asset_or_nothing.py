"""Two-asset asset-or-nothing digital (two_asset_asset_or_nothing)."""

import math
import random

import pytest

from quantforge import (
    two_asset_asset_or_nothing,
    asset_or_nothing,
    cash_or_nothing,
    OptionType,
)


S1, S2, K1, K2, T, R = 100.0, 95.0, 105.0, 90.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


def _mc(cond1, cond2, rho, seed, n=800_000):
    rng = random.Random(seed)
    d1 = (R - 0.5 * SIG1 * SIG1) * T
    d2 = (R - 0.5 * SIG2 * SIG2) * T
    v1 = SIG1 * math.sqrt(T)
    v2 = SIG2 * math.sqrt(T)
    c2 = math.sqrt(1.0 - rho * rho)
    disc = math.exp(-R * T)
    tot = 0.0
    for _ in range(n):
        z1 = rng.gauss(0.0, 1.0)
        z2 = rng.gauss(0.0, 1.0)
        a = S1 * math.exp(d1 + v1 * z1)
        b = S2 * math.exp(d2 + v2 * (rho * z1 + c2 * z2))
        o1 = (a > K1) if cond1 == "above" else (a < K1)
        o2 = (b > K2) if cond2 == "above" else (b < K2)
        if o1 and o2:
            tot += a
    return disc * tot / n


@pytest.mark.slow
@pytest.mark.parametrize("c1,c2", [("above", "above"), ("above", "below"),
                                   ("below", "above"), ("below", "below")])
def test_quadrants_match_mc(c1, c2):
    closed = two_asset_asset_or_nothing(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                        c1, c2)
    mc = _mc(c1, c2, RHO, seed=hash((c1, c2)) % 10_000, n=1_500_000)
    # Asset-weighted payoff -> higher MC variance; a ~0.15 band covers 3 SE.
    assert closed == pytest.approx(mc, abs=0.15)


def test_quadrants_sum_to_forward():
    # Asset 1 is always delivered on some quadrant -> sum = S1 e^{-q1 t}.
    total = sum(two_asset_asset_or_nothing(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                           c1, c2)
                for c1 in ("above", "below") for c2 in ("above", "below"))
    assert total == pytest.approx(S1, abs=1e-8)


def test_zero_correlation_factorizes():
    # rho = 0: pays S1 iff S1>K1 and S2>K2 = (asset-or-nothing on S1) * P(S2>K2).
    disc = math.exp(-R * T)
    joint = two_asset_asset_or_nothing(S1, S2, K1, K2, T, R, SIG1, SIG2, 0.0)
    a1 = asset_or_nothing(S1, K1, T, R, SIG1, OptionType.CALL)
    p2 = cash_or_nothing(S2, K2, T, R, SIG2, OptionType.CALL) / disc
    assert joint == pytest.approx(a1 * p2, abs=1e-8)


def test_dividend_forward_sum():
    q1 = 0.03
    total = sum(two_asset_asset_or_nothing(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                           c1, c2, q1=q1)
                for c1 in ("above", "below") for c2 in ("above", "below"))
    assert total == pytest.approx(S1 * math.exp(-q1 * T), abs=1e-8)


def test_bad_condition_raises():
    with pytest.raises(ValueError):
        two_asset_asset_or_nothing(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO,
                                   cond2="maybe")


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        two_asset_asset_or_nothing(S1, S2, K1, K2, T, R, SIG1, SIG2, 2.0)
