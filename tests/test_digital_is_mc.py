"""Cash-or-nothing digital by importance sampling (digital_is_mc)."""

import math
import random

import pytest

from quantforge import (
    OptionType,
    digital_is_mc,
    cash_or_nothing,
)


S, T, R, SIG = 100.0, 1.0, 0.03, 0.2


def _plain_indicator_se(K, n=40_000, seed=0):
    rng = random.Random(seed)
    disc = math.exp(-R * T)
    drift = (R - 0.5 * SIG * SIG) * T
    vol = SIG * math.sqrt(T)
    xs = []
    for _ in range(n):
        z = rng.gauss(0.0, 1.0)
        sT = S * math.exp(drift + vol * z)
        xs.append(disc if sT > K else 0.0)
    m = sum(xs) / n
    var = sum((x - m) ** 2 for x in xs) / (n - 1)
    return math.sqrt(var / n)


def test_deep_otm_call_matches_closed_form():
    K = 160.0
    cf = cash_or_nothing(S, K, T, R, SIG, OptionType.CALL)
    mc = digital_is_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=60_000, seed=1)
    assert mc.price == pytest.approx(cf, abs=3.0 * mc.std_error)


def test_deep_otm_put_matches_closed_form():
    K = 55.0
    cf = cash_or_nothing(S, K, T, R, SIG, OptionType.PUT)
    mc = digital_is_mc(S, K, T, R, SIG, OptionType.PUT, n_paths=60_000, seed=2)
    assert mc.price == pytest.approx(cf, abs=3.0 * mc.std_error)


def test_importance_sampling_beats_plain_indicator():
    K = 160.0
    mc = digital_is_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=40_000, seed=3)
    plain_se = _plain_indicator_se(K, n=40_000, seed=3)
    assert mc.std_error < 0.5 * plain_se


def test_atm_still_correct():
    cf = cash_or_nothing(S, 100.0, T, R, SIG, OptionType.CALL)
    mc = digital_is_mc(S, 100.0, T, R, SIG, OptionType.CALL, n_paths=60_000, seed=4)
    assert mc.price == pytest.approx(cf, abs=3.0 * mc.std_error)


def test_cash_scales_linearly():
    K = 150.0
    mc = digital_is_mc(S, K, T, R, SIG, OptionType.CALL, cash=10.0,
                       n_paths=60_000, seed=5)
    cf = cash_or_nothing(S, K, T, R, SIG, OptionType.CALL, cash=10.0)
    assert mc.price == pytest.approx(cf, abs=3.0 * mc.std_error)


def test_dividend_carry():
    K = 150.0
    b = R - 0.04
    cf = cash_or_nothing(S, K, T, R, SIG, OptionType.CALL, b=b)
    mc = digital_is_mc(S, K, T, R, SIG, OptionType.CALL, b=b,
                       n_paths=60_000, seed=6)
    assert mc.price == pytest.approx(cf, abs=3.0 * mc.std_error)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        digital_is_mc(-1, 160.0, T, R, SIG)
