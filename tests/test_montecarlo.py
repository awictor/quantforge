"""Tests for the Monte Carlo engine.

Seeds are fixed so runs are deterministic. We check convergence to closed-form
values and, critically, that the control-variate estimator has a much smaller
standard error than plain Monte Carlo.
"""

import math

import pytest

from quantforge import (
    call_price, put_price, geometric_asian, OptionType,
    european_mc, arithmetic_asian_mc,
)


# --- European MC converges to BSM ---
@pytest.mark.parametrize("ot,ref", [
    (OptionType.CALL, None),
    (OptionType.PUT, None),
])
def test_european_mc_matches_bsm(ot, ref):
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.2
    res = european_mc(S, K, t, r, sigma, ot, n_paths=200_000, seed=42)
    exact = call_price(S, K, t, r, sigma) if ot is OptionType.CALL else put_price(S, K, t, r, sigma)
    # Within ~3 standard errors of the analytic price.
    assert abs(res.price - exact) < 3 * res.std_error + 1e-9
    # And absolutely close given the path count.
    assert res.price == pytest.approx(exact, abs=0.05)


def test_mc_result_confidence_interval():
    res = european_mc(100, 100, 1.0, 0.05, 0.2, OptionType.CALL,
                      n_paths=50_000, seed=7)
    lo, hi = res.confidence_interval()
    assert lo < res.price < hi
    assert hi - lo == pytest.approx(2 * 1.96 * res.std_error)


def test_antithetic_reduces_error():
    kw = dict(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
              option_type=OptionType.CALL, n_paths=40_000, seed=123)
    plain = european_mc(antithetic=False, **kw)
    anti = european_mc(antithetic=True, **kw)
    # Antithetic sampling should not increase the error; usually cuts it.
    assert anti.std_error <= plain.std_error * 1.05


# --- Arithmetic Asian: control variate ---
def test_control_variate_shrinks_std_error():
    kw = dict(S=100, K=100, t=1.0, r=0.05, sigma=0.3,
              option_type=OptionType.CALL, n_steps=50, n_paths=20_000, seed=99)
    plain = arithmetic_asian_mc(control_variate=False, **kw)
    cv = arithmetic_asian_mc(control_variate=True, **kw)
    # The geometric control is almost perfectly correlated -> big SE drop.
    assert cv.std_error < plain.std_error / 5
    # Both estimate the same price, within plain MC's error band.
    assert abs(cv.price - plain.price) < 3 * plain.std_error


def test_arithmetic_asian_above_geometric():
    # By the AM-GM inequality the arithmetic-average call is worth at least the
    # geometric-average call.
    S, K, t, r, sigma = 100, 100, 1.0, 0.05, 0.3
    cv = arithmetic_asian_mc(S, K, t, r, sigma, OptionType.CALL,
                             n_steps=50, n_paths=40_000, seed=5)
    geo = geometric_asian(S, K, t, r, sigma, OptionType.CALL)
    assert cv.price > geo - 3 * cv.std_error
    assert cv.price >= geo - 1e-3  # arithmetic >= geometric (up to MC noise)


def test_seed_is_reproducible():
    a = arithmetic_asian_mc(100, 100, 1.0, 0.05, 0.3, n_steps=20,
                            n_paths=5_000, seed=2024)
    b = arithmetic_asian_mc(100, 100, 1.0, 0.05, 0.3, n_steps=20,
                            n_paths=5_000, seed=2024)
    assert a.price == b.price
    assert a.std_error == b.std_error
