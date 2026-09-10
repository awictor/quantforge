"""Tests for the trinomial lattice and Richardson extrapolation."""

import pytest

from quantforge import (
    trinomial_price, richardson_american, american_price,
    call_price, put_price, OptionType,
)


def test_european_trinomial_converges_to_bsm_call():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    tri = trinomial_price(S, K, t, r, sigma, OptionType.CALL, steps=400, american=False)
    assert tri == pytest.approx(call_price(S, K, t, r, sigma), abs=1e-2)


def test_european_trinomial_converges_to_bsm_put():
    S, K, t, r, sigma = 100, 105, 0.5, 0.03, 0.3
    tri = trinomial_price(S, K, t, r, sigma, OptionType.PUT, steps=400, american=False)
    assert tri == pytest.approx(put_price(S, K, t, r, sigma), abs=1e-2)


@pytest.mark.slow
def test_american_trinomial_matches_binomial():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    tri = trinomial_price(S, K, t, r, sigma, OptionType.PUT, steps=400, american=True)
    binom = american_price(S, K, t, r, sigma, OptionType.PUT, steps=2000)
    assert tri == pytest.approx(binom, abs=1e-2)


def test_american_call_no_dividend_equals_european():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    tri = trinomial_price(S, K, t, r, sigma, OptionType.CALL, steps=400, american=True)
    assert tri == pytest.approx(call_price(S, K, t, r, sigma), abs=1e-2)


def test_american_put_early_exercise_premium():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    am = trinomial_price(S, K, t, r, sigma, OptionType.PUT, steps=400, american=True)
    eu = put_price(S, K, t, r, sigma)
    assert am > eu  # early exercise is worth something here


@pytest.mark.slow
def test_richardson_more_accurate_than_plain_trinomial():
    # Reference: a very fine binomial tree.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    ref = american_price(S, K, t, r, sigma, OptionType.PUT, steps=5000)
    plain = trinomial_price(S, K, t, r, sigma, OptionType.PUT, steps=50, american=True)
    rich = richardson_american(S, K, t, r, sigma, OptionType.PUT, steps=50)
    assert abs(rich - ref) < abs(plain - ref)


@pytest.mark.slow
def test_dividend_put_matches_binomial():
    # Carry b = r - q with a dividend yield.
    S, K, t, r, sigma, b = 100, 100, 1.0, 0.06, 0.25, 0.06 - 0.03
    tri = trinomial_price(S, K, t, r, sigma, OptionType.PUT, b=b, steps=400)
    binom = american_price(S, K, t, r, sigma, OptionType.PUT, b=b, steps=2000)
    assert tri == pytest.approx(binom, abs=2e-2)


def test_zero_time_is_intrinsic():
    assert trinomial_price(120, 100, 0.0, 0.05, 0.2, OptionType.CALL) == pytest.approx(20.0)
    assert trinomial_price(80, 100, 0.0, 0.05, 0.2, OptionType.PUT) == pytest.approx(20.0)
