"""Leisen-Reimer binomial tree: fast convergence and American consistency."""

import pytest

from quantforge import (
    OptionType,
    leisen_reimer_price as lr,
    leisen_reimer_greeks,
    american_price as crr,
    bjerksund_stensland as bs,
    call_price,
    put_price,
)


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.2


def test_european_matches_black_scholes_at_low_steps():
    # Even ~51 steps lands within a fraction of a cent of Black-Scholes.
    assert lr(S, K, T, R, SIGMA, OptionType.CALL, steps=51) == pytest.approx(
        call_price(S, K, T, R, SIGMA), abs=2e-3)
    assert lr(S, K, T, R, SIGMA, OptionType.PUT, steps=51) == pytest.approx(
        put_price(S, K, T, R, SIGMA), abs=2e-3)


def test_converges_faster_than_crr():
    # At the same modest step count LR is far closer to Black-Scholes than CRR.
    bp = put_price(S, K, T, R, SIGMA)
    n = 51
    lr_err = abs(lr(S, K, T, R, SIGMA, OptionType.PUT, steps=n) - bp)
    crr_err = abs(crr(S, K, T, R, SIGMA, OptionType.PUT, steps=n) - bp)
    assert lr_err < crr_err
    assert lr_err < 1e-3


def test_error_shrinks_with_steps():
    bp = put_price(S, K, T, R, SIGMA)
    e_small = abs(lr(S, K, T, R, SIGMA, OptionType.PUT, steps=21) - bp)
    e_large = abs(lr(S, K, T, R, SIGMA, OptionType.PUT, steps=201) - bp)
    assert e_large < e_small


def test_even_steps_forced_odd():
    # An even request is bumped to odd; 50 and 51 give the same price.
    a = lr(S, K, T, R, SIGMA, OptionType.CALL, steps=50)
    b = lr(S, K, T, R, SIGMA, OptionType.CALL, steps=51)
    assert a == pytest.approx(b, abs=1e-12)


AMERICAN_CASES = [
    (42, 40, 0.75, 0.04, 0.35, OptionType.CALL, -0.04),
    (100, 100, 0.5, 0.05, 0.30, OptionType.PUT, 0.05),
    (110, 100, 0.5, 0.10, 0.25, OptionType.PUT, 0.10),
    (90, 100, 1.0, 0.08, 0.20, OptionType.PUT, 0.08),
]


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,sigma,ot,b", AMERICAN_CASES)
def test_american_matches_converged_crr_with_dividends(S, K, t, r, sigma, ot, b):
    # American convergence is slower than the European O(1/n^2) (the smooth-
    # payoff assumption behind Peizer-Pratt breaks at the exercise boundary), so
    # use a larger step count; LR still reaches the 3000-step CRR tree to ~1 cent.
    lr_val = lr(S, K, t, r, sigma, ot, b=b, steps=601, american=True)
    tree = crr(S, K, t, r, sigma, ot, b=b, steps=3000)
    assert lr_val == pytest.approx(tree, abs=1.5e-2)


def test_no_dividend_american_call_equals_european():
    # b = r: never early-exercise a call, so the American LR equals European.
    am = lr(S, K, T, R, SIGMA, OptionType.CALL, steps=101, american=True)
    assert am == pytest.approx(call_price(S, K, T, R, SIGMA), abs=2e-3)


def test_american_put_above_european():
    am = lr(S, K, T, R, SIGMA, OptionType.PUT, steps=101, american=True)
    eu = lr(S, K, T, R, SIGMA, OptionType.PUT, steps=101, american=False)
    assert am > eu


def test_greeks_delta_gamma_reasonable():
    g = leisen_reimer_greeks(S, K, T, R, SIGMA, OptionType.CALL, steps=101)
    assert 0.0 < g["delta"] < 1.0
    assert g["gamma"] > 0.0
    assert g["vega"] > 0.0


def test_zero_time_is_intrinsic():
    assert lr(120, 100, 0.0, 0.05, 0.2, OptionType.CALL) == pytest.approx(20.0)
    assert lr(80, 100, 0.0, 0.05, 0.2, OptionType.PUT) == pytest.approx(20.0)
