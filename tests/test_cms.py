"""CMS convexity adjustment: standard model vs static replication."""

import math

import pytest

from quantforge import (
    cms_adjustment_standard,
    cms_rate_convexity_replication,
    cms_rate,
)


F, SIG, EXP, TENOR = 0.04, 0.2, 5.0, 10.0


def test_zero_vol_no_adjustment():
    assert cms_adjustment_standard(F, 0.0, EXP, TENOR) == 0.0


def test_adjustment_positive():
    assert cms_adjustment_standard(F, SIG, EXP, TENOR) > 0.0


def test_replication_flat_smile_matches_standard():
    std = cms_adjustment_standard(F, SIG, EXP, TENOR)
    rep = cms_rate_convexity_replication(F, EXP, TENOR, lambda K: SIG)
    assert rep == pytest.approx(std, abs=1e-5)


def test_replication_matches_linear_tsr_monte_carlo():
    # E_pay[S_T] - S0 = G * Var_A(S_T) under the linear TSR; MC-check the value.
    import random
    from quantforge.cms import _level_G
    G = _level_G(F, TENOR, 1.0, 0.0)
    rng = random.Random(1)
    n = 500_000
    num = den = 0.0
    for _ in range(n):
        z = rng.gauss(0, 1)
        ST = F * math.exp(-0.5 * SIG * SIG * EXP + SIG * math.sqrt(EXP) * z)
        al = 1.0 + G * (ST - F)
        num += ST * al
        den += al
    mc_ca = num / den - F
    assert cms_adjustment_standard(F, SIG, EXP, TENOR) == pytest.approx(mc_ca,
                                                                        abs=2e-4)


def test_smile_raises_adjustment():
    def smile(K):
        return SIG + 0.1 * (K - F) ** 2 / (F * F)   # convex smile
    flat = cms_rate_convexity_replication(F, EXP, TENOR, lambda K: SIG)
    smiled = cms_rate_convexity_replication(F, EXP, TENOR, smile)
    assert smiled > flat


def test_cms_rate_adds_adjustment():
    r = cms_rate(F, SIG, EXP, TENOR)
    assert r == pytest.approx(F + cms_adjustment_standard(F, SIG, EXP, TENOR))
    assert r > F


def test_pay_lag_increases_adjustment():
    natural = cms_adjustment_standard(F, SIG, EXP, TENOR, pay_lag=0.0)
    lagged = cms_adjustment_standard(F, SIG, EXP, TENOR, pay_lag=0.5)
    assert lagged > natural
