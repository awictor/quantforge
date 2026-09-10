"""Tests for Longstaff-Schwartz least-squares Monte Carlo."""

import pytest

from quantforge import bermudan_lsm, american_price, put_price, call_price, OptionType
from quantforge.lsm import _solve_normal_equations


def test_least_squares_fits_a_line():
    # y = 2x + 1 sampled exactly -> beta = [1, 2].
    xs = [0.0, 1.0, 2.0, 3.0, 4.0]
    X = [[1.0, x] for x in xs]
    y = [1 + 2 * x for x in xs]
    beta = _solve_normal_equations(X, y)
    assert beta[0] == pytest.approx(1.0, abs=1e-9)
    assert beta[1] == pytest.approx(2.0, abs=1e-9)


@pytest.mark.slow
def test_bermudan_put_matches_binomial_american():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    lsm = bermudan_lsm(S, K, t, r, sigma, OptionType.PUT, n_steps=50,
                       n_paths=40_000, seed=1)
    tree = american_price(S, K, t, r, sigma, OptionType.PUT, steps=2000)
    # LSM with 50 dates approaches the American value from slightly below.
    assert lsm == pytest.approx(tree, abs=0.15)


@pytest.mark.slow
def test_bermudan_put_at_least_european():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    lsm = bermudan_lsm(S, K, t, r, sigma, OptionType.PUT, n_steps=50,
                       n_paths=40_000, seed=3)
    eu = put_price(S, K, t, r, sigma)
    # Early exercise has value, but the LSM estimate is mildly biased low, so
    # allow a small tolerance below the European price.
    assert lsm >= eu - 0.05


@pytest.mark.slow
def test_no_dividend_call_close_to_european():
    # American call without dividends == European call; LSM should be near it.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    lsm = bermudan_lsm(S, K, t, r, sigma, OptionType.CALL, n_steps=50,
                       n_paths=40_000, seed=5)
    eu = call_price(S, K, t, r, sigma)
    assert lsm == pytest.approx(eu, abs=0.2)


def test_reproducible_with_seed():
    a = bermudan_lsm(100, 100, 0.5, 0.05, 0.25, OptionType.PUT, n_steps=20,
                     n_paths=5_000, seed=7)
    b = bermudan_lsm(100, 100, 0.5, 0.05, 0.25, OptionType.PUT, n_steps=20,
                     n_paths=5_000, seed=7)
    assert a == b


def test_rejects_bad_steps():
    with pytest.raises(ValueError):
        bermudan_lsm(100, 100, 1.0, 0.05, 0.2, n_steps=0)
