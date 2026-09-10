"""American option pricing via Kim's (1990) integral equation vs binomial."""

import pytest

from quantforge import (
    OptionType,
    kim_american_put,
    kim_american_call,
    kim_exercise_boundary,
    american_price as crr,
    put_price,
    call_price,
)


PUT_CASES = [
    (100, 100, 1.0, 0.05, 0.20, 0.0),
    (90, 100, 0.5, 0.05, 0.30, 0.0),
    (100, 100, 1.0, 0.05, 0.20, 0.03),
    (110, 100, 1.0, 0.08, 0.25, 0.04),
]


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,sigma,q", PUT_CASES)
def test_put_matches_binomial(S, K, t, r, sigma, q):
    kim = kim_american_put(S, K, t, r, sigma, q=q, n_steps=200)
    tree = crr(S, K, t, r, sigma, OptionType.PUT, b=r - q, steps=4000)
    assert kim == pytest.approx(tree, abs=3e-3)


@pytest.mark.slow
@pytest.mark.parametrize("S,K,t,r,sigma,q", [
    (100, 100, 1.0, 0.05, 0.20, 0.06),
    (110, 100, 0.5, 0.03, 0.25, 0.05),
])
def test_dividend_call_matches_binomial(S, K, t, r, sigma, q):
    kim = kim_american_call(S, K, t, r, sigma, q=q, n_steps=200)
    tree = crr(S, K, t, r, sigma, OptionType.CALL, b=r - q, steps=4000)
    assert kim == pytest.approx(tree, abs=3e-3)


def test_no_dividend_call_equals_european():
    ac = kim_american_call(100, 100, 1.0, 0.05, 0.2, q=0.0)
    assert ac == pytest.approx(call_price(100, 100, 1.0, 0.05, 0.2), abs=1e-9)


def test_put_at_least_european_and_intrinsic():
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    am = kim_american_put(S, K, t, r, sigma, n_steps=120)
    eu = put_price(S, K, t, r, sigma)
    assert am >= eu - 1e-6
    assert am > eu  # early exercise has value with no dividends on a put


def test_deep_itm_put_is_intrinsic():
    # Below the current exercise boundary the option is exercised now.
    v = kim_american_put(50, 100, 1.0, 0.05, 0.2, n_steps=80)
    assert v == pytest.approx(50.0, abs=1e-6)


def test_boundary_increasing_to_strike():
    B = kim_exercise_boundary(100, 1.0, 0.05, 0.2, q=0.0, n_steps=60)
    assert B[-1] == pytest.approx(100.0)          # B(T) = K for q = 0
    assert all(B[i] <= B[i + 1] + 1e-9 for i in range(len(B) - 1))
    assert B[0] < 100.0                            # exercise below strike early


def test_zero_time_is_intrinsic():
    assert kim_american_put(80, 100, 0.0, 0.05, 0.2) == pytest.approx(20.0)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        kim_american_put(-1, 100, 1.0, 0.05, 0.2)
    with pytest.raises(ValueError):
        kim_american_put(100, 100, 1.0, 0.05, 0.0)  # sigma <= 0
