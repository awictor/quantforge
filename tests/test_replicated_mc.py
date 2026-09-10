"""Batched-replication honest standard error (replicated_mc)."""

import math

import pytest

from quantforge import (
    OptionType,
    replicated_mc,
    european_mc,
    european_stratified_mc,
    call_price,
)


S, K, T, R, SIG = 100.0, 105.0, 1.0, 0.04, 0.25


def test_price_unbiased_for_stratified():
    bs = call_price(S, K, T, R, SIG)
    rep = replicated_mc(
        lambda s: european_stratified_mc(S, K, T, R, SIG, OptionType.CALL,
                                          n_strata=200, n_per=10, seed=s),
        n_batches=40)
    assert rep.price == pytest.approx(bs, abs=3.0 * rep.std_error)


def test_n_paths_is_batch_count():
    rep = replicated_mc(
        lambda s: european_stratified_mc(S, K, T, R, SIG, seed=s),
        n_batches=25)
    assert rep.n_paths == 25


def test_accepts_float_estimator():
    # An estimator returning a bare float (not an MCResult) also works.
    rep = replicated_mc(lambda s: european_mc(S, K, T, R, SIG, n_paths=2000,
                                              seed=s).price,
                        n_batches=20)
    assert rep.price > 0.0
    assert rep.std_error > 0.0


@pytest.mark.slow
def test_honest_se_calibrated_for_plain_estimator():
    # For a plain i.i.d. estimator, the across-batch SE of R batches of n paths
    # each must match a single-run i.i.d. SE / sqrt(R) up to sampling noise.
    R_BATCH, N = 60, 5000
    rep = replicated_mc(
        lambda s: european_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=N, seed=s),
        n_batches=R_BATCH)
    one = european_mc(S, K, T, R, SIG, OptionType.CALL, n_paths=N, seed=999)
    expected = one.std_error / math.sqrt(R_BATCH)
    assert rep.std_error == pytest.approx(expected, rel=0.4)


@pytest.mark.slow
def test_honest_se_reveals_lhs_reduction():
    # The LHS spread estimator reports a naive i.i.d. std_error that cannot see
    # its variance reduction; the honest (across-seed) SE is far smaller than a
    # plain two-asset MC's honest SE at equal paths.
    from quantforge import spread_option_lhs_mc

    S1, S2, KS = 100.0, 95.0, 5.0
    sig1, sig2, rho = 0.25, 0.30, 0.5

    def plain_spread(seed, n=4000):
        import random
        rng = random.Random(seed)
        d1 = (R - 0.5 * sig1 * sig1) * T
        d2 = (R - 0.5 * sig2 * sig2) * T
        v1 = sig1 * math.sqrt(T)
        v2 = sig2 * math.sqrt(T)
        c2 = math.sqrt(1.0 - rho * rho)
        disc = math.exp(-R * T)
        tot = 0.0
        for _ in range(n):
            z1 = rng.gauss(0.0, 1.0)
            z2 = rng.gauss(0.0, 1.0)
            a = S1 * math.exp(d1 + v1 * z1)
            bb = S2 * math.exp(d2 + v2 * (rho * z1 + c2 * z2))
            tot += disc * max(a - bb - KS, 0.0)
        return tot / n

    lhs = replicated_mc(
        lambda s: spread_option_lhs_mc(S1, S2, KS, T, R, sig1, sig2, rho,
                                       n_paths=4000, seed=s),
        n_batches=40)
    plain = replicated_mc(plain_spread, n_batches=40)
    assert lhs.std_error < 0.6 * plain.std_error


def test_bad_n_batches_raises():
    with pytest.raises(ValueError):
        replicated_mc(lambda s: european_mc(S, K, T, R, SIG, seed=s), n_batches=1)
