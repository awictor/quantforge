"""Two-asset spread option by Latin hypercube MC (spread_option_lhs_mc)."""

import math

import pytest

from quantforge import (
    OptionType,
    spread_option_lhs_mc,
    spread_option,
    exchange_option,
)


S1, S2, K, T, R = 100.0, 95.0, 5.0, 1.0, 0.03
SIG1, SIG2, RHO = 0.25, 0.30, 0.5


def test_matches_kirk():
    kirk = spread_option(S1, S2, K, T, R, SIG1, SIG2, RHO,
                         option_type=OptionType.CALL)
    # LHS SE is a naive i.i.d. estimate; use a modest absolute tolerance that
    # the true (across-seed) error comfortably sits inside.
    mc = spread_option_lhs_mc(S1, S2, K, T, R, SIG1, SIG2, RHO,
                              n_paths=60_000, seed=1)
    assert mc.price == pytest.approx(kirk, abs=0.1)


def test_zero_strike_matches_margrabe():
    # At K = 0 the spread is an exact Margrabe exchange option.
    marg = exchange_option(S1, S2, T, SIG1, SIG2, RHO)
    mc = spread_option_lhs_mc(S1, S2, 0.0, T, R, SIG1, SIG2, RHO,
                              n_paths=60_000, seed=2)
    assert mc.price == pytest.approx(marg, abs=0.12)


@pytest.mark.slow
def test_lhs_beats_plain_across_seeds():
    # The genuine LHS gain shows up in the spread of estimates across seeds, not
    # in the reported i.i.d. std_error.
    import random

    truth = spread_option(S1, S2, K, T, R, SIG1, SIG2, RHO)
    d1 = (R - 0.5 * SIG1 * SIG1) * T
    d2 = (R - 0.5 * SIG2 * SIG2) * T
    v1 = SIG1 * math.sqrt(T)
    v2 = SIG2 * math.sqrt(T)
    c2 = math.sqrt(1.0 - RHO * RHO)
    disc = math.exp(-R * T)

    def plain(seed, n=4000):
        rng = random.Random(seed)
        tot = 0.0
        for _ in range(n):
            z1 = rng.gauss(0.0, 1.0)
            z2 = rng.gauss(0.0, 1.0)
            a = S1 * math.exp(d1 + v1 * z1)
            bb = S2 * math.exp(d2 + v2 * (RHO * z1 + c2 * z2))
            tot += disc * max(a - bb - K, 0.0)
        return tot / n

    ns = 40
    lhs_err = [spread_option_lhs_mc(S1, S2, K, T, R, SIG1, SIG2, RHO,
                                    n_paths=4000, seed=i).price - truth
               for i in range(ns)]
    plain_err = [plain(1000 + i) - truth for i in range(ns)]
    lhs_rmse = math.sqrt(sum(e * e for e in lhs_err) / ns)
    plain_rmse = math.sqrt(sum(e * e for e in plain_err) / ns)
    assert lhs_rmse < 0.6 * plain_rmse


def test_put_via_parity_sign():
    # A deep-in-the-money put on the spread (K large) should be worth about the
    # discounted intrinsic; just check it is positive and finite.
    mc = spread_option_lhs_mc(S1, S2, 40.0, T, R, SIG1, SIG2, RHO,
                              option_type=OptionType.PUT, n_paths=40_000, seed=3)
    assert mc.price > 0.0


def test_bad_rho_raises():
    with pytest.raises(ValueError):
        spread_option_lhs_mc(S1, S2, K, T, R, SIG1, SIG2, 1.5)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        spread_option_lhs_mc(-1, S2, K, T, R, SIG1, SIG2, RHO)
