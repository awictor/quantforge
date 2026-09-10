"""Two-asset correlated cash-or-nothing digital (two_asset_digital)."""

import math
import random

import pytest

from quantforge import two_asset_digital, cash_or_nothing, OptionType


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
    hit = 0
    for _ in range(n):
        z1 = rng.gauss(0.0, 1.0)
        z2 = rng.gauss(0.0, 1.0)
        a = S1 * math.exp(d1 + v1 * z1)
        b = S2 * math.exp(d2 + v2 * (rho * z1 + c2 * z2))
        o1 = (a > K1) if cond1 == "above" else (a < K1)
        o2 = (b > K2) if cond2 == "above" else (b < K2)
        if o1 and o2:
            hit += 1
    return disc * hit / n


@pytest.mark.slow
@pytest.mark.parametrize("c1,c2", [("above", "above"), ("above", "below"),
                                   ("below", "above"), ("below", "below")])
def test_quadrants_match_mc(c1, c2):
    closed = two_asset_digital(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO, c1, c2)
    mc = _mc(c1, c2, RHO, seed=hash((c1, c2)) % 10_000)
    assert closed == pytest.approx(mc, abs=3e-3)


def test_quadrants_sum_to_discount():
    # The four exhaustive conditions partition all outcomes -> sum = e^{-rt}.
    total = sum(two_asset_digital(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO, c1, c2)
                for c1 in ("above", "below") for c2 in ("above", "below"))
    assert total == pytest.approx(math.exp(-R * T), abs=1e-9)


def test_zero_correlation_factorizes():
    # At rho = 0 the joint digital is the product of the two single-asset
    # digitals (divided by the extra discount factor each carries).
    disc = math.exp(-R * T)
    joint = two_asset_digital(S1, S2, K1, K2, T, R, SIG1, SIG2, 0.0)
    d1 = cash_or_nothing(S1, K1, T, R, SIG1, OptionType.CALL)
    d2 = cash_or_nothing(S2, K2, T, R, SIG2, OptionType.CALL)
    assert joint == pytest.approx(d1 * d2 / disc, abs=1e-9)


def test_cash_scales_linearly():
    a = two_asset_digital(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO, cash=1.0)
    b = two_asset_digital(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO, cash=7.5)
    assert b == pytest.approx(7.5 * a, abs=1e-12)


def test_bad_condition_raises():
    with pytest.raises(ValueError):
        two_asset_digital(S1, S2, K1, K2, T, R, SIG1, SIG2, RHO, cond1="sideways")


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        two_asset_digital(S1, S2, K1, K2, T, R, SIG1, SIG2, 1.5)
