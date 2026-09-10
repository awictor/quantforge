"""Two-asset gap option gated by a second asset (two_asset_gap_option)."""

import math
import random

import pytest

from quantforge import (
    two_asset_gap_option,
    correlation_option,
    OptionType,
)


S1, S2, KTRIG, KPAY, K2, T, R = 100.0, 95.0, 100.0, 105.0, 90.0, 1.0, 0.05
SIG1, SIG2, RHO = 0.2, 0.3, 0.4


def _mc(ot, cond2, seed, ktrig=KTRIG, kpay=KPAY, n=800_000):
    rng = random.Random(seed)
    d1 = (R - 0.5 * SIG1 * SIG1) * T
    d2 = (R - 0.5 * SIG2 * SIG2) * T
    v1 = SIG1 * math.sqrt(T)
    v2 = SIG2 * math.sqrt(T)
    c2 = math.sqrt(1.0 - RHO * RHO)
    disc = math.exp(-R * T)
    tot = 0.0
    for _ in range(n):
        z1 = rng.gauss(0.0, 1.0)
        z2 = rng.gauss(0.0, 1.0)
        a = S1 * math.exp(d1 + v1 * z1)
        b = S2 * math.exp(d2 + v2 * (RHO * z1 + c2 * z2))
        gate = (b > K2) if cond2 == "above" else (b < K2)
        if not gate:
            continue
        if ot is OptionType.CALL:
            if a > ktrig:
                tot += a - kpay
        else:
            if a < ktrig:
                tot += kpay - a
    return disc * tot / n


@pytest.mark.slow
def test_call_matches_mc():
    closed = two_asset_gap_option(S1, S2, KTRIG, KPAY, K2, T, R, SIG1, SIG2, RHO,
                                  OptionType.CALL, "above")
    assert closed == pytest.approx(_mc(OptionType.CALL, "above", 1), abs=0.02)


@pytest.mark.slow
def test_put_matches_mc():
    closed = two_asset_gap_option(S1, S2, KTRIG, 95.0, K2, T, R, SIG1, SIG2, RHO,
                                  OptionType.PUT, "above")
    assert closed == pytest.approx(
        _mc(OptionType.PUT, "above", 2, kpay=95.0), abs=0.02)


def test_reduces_to_correlation_option():
    # Equal trigger and payoff strikes -> plain correlation option.
    g = two_asset_gap_option(S1, S2, 100.0, 100.0, K2, T, R, SIG1, SIG2, RHO,
                             OptionType.CALL, "above")
    c = correlation_option(S1, S2, 100.0, K2, T, R, SIG1, SIG2, RHO,
                           OptionType.CALL, "above")
    assert g == pytest.approx(c, abs=1e-12)


def test_higher_payoff_strike_lowers_call():
    # A larger payoff strike subtracts more per unit -> lower call value.
    lo = two_asset_gap_option(S1, S2, KTRIG, 100.0, K2, T, R, SIG1, SIG2, RHO,
                              OptionType.CALL, "above")
    hi = two_asset_gap_option(S1, S2, KTRIG, 110.0, K2, T, R, SIG1, SIG2, RHO,
                              OptionType.CALL, "above")
    assert hi < lo


def test_bad_cond_raises():
    with pytest.raises(ValueError):
        two_asset_gap_option(S1, S2, KTRIG, KPAY, K2, T, R, SIG1, SIG2, RHO,
                             cond2="maybe")


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        two_asset_gap_option(S1, S2, KTRIG, KPAY, K2, T, R, SIG1, SIG2, 1.5)
