"""Two-asset basket option by Latin hypercube MC (basket_option_lhs_mc)."""

import pytest

from quantforge import (
    OptionType,
    basket_option_lhs_mc,
    basket_option,
    replicated_mc,
)


SPOTS, W = (100.0, 90.0), (0.6, 0.4)
K, T, R = 95.0, 1.0, 0.03
SIG, CORR = (0.2, 0.3), 0.4


def test_call_near_moment_match():
    # basket_option is a Levy moment-match (approximate); LHS is the unbiased
    # MC. They agree to a few cents for typical vols.
    mm = basket_option(SPOTS, W, K, T, R, SIG, CORR, option_type=OptionType.CALL)
    mc = basket_option_lhs_mc(SPOTS, W, K, T, R, SIG, CORR,
                              n_paths=60_000, seed=1)
    assert mc.price == pytest.approx(mm, abs=0.1)


def test_put_near_moment_match():
    mm = basket_option(SPOTS, W, K, T, R, SIG, CORR, option_type=OptionType.PUT)
    mc = basket_option_lhs_mc(SPOTS, W, K, T, R, SIG, CORR,
                              option_type=OptionType.PUT, n_paths=60_000, seed=2)
    assert mc.price == pytest.approx(mm, abs=0.1)


def test_honest_se_via_replication():
    # Wrapped in replicated_mc the across-seed SE is honest; the moment-match
    # price sits within a few of those SEs.
    mm = basket_option(SPOTS, W, K, T, R, SIG, CORR, option_type=OptionType.CALL)
    rep = replicated_mc(
        lambda s: basket_option_lhs_mc(SPOTS, W, K, T, R, SIG, CORR,
                                       n_paths=8000, seed=s),
        n_batches=40)
    assert rep.price == pytest.approx(mm, abs=5.0 * rep.std_error)


def test_single_asset_weight_reduces_to_bs():
    # Weight all mass on asset 1 (w2 = 0): a plain Black-Scholes call on 100.
    from quantforge import call_price
    bs = call_price(100.0, K, T, R, 0.2)
    mc = basket_option_lhs_mc((100.0, 90.0), (1.0, 0.0), K, T, R, (0.2, 0.3),
                              CORR, n_paths=60_000, seed=3)
    assert mc.price == pytest.approx(bs, abs=0.05)


def test_bad_corr_raises():
    with pytest.raises(ValueError):
        basket_option_lhs_mc(SPOTS, W, K, T, R, SIG, 1.5)


def test_bad_dimension_raises():
    with pytest.raises(ValueError):
        basket_option_lhs_mc((100.0,), (1.0,), K, T, R, (0.2,), CORR)
