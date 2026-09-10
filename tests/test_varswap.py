"""Tests for variance/volatility swap fair-strike replication.

Under Black-Scholes (flat vol sigma) the fair variance-swap strike equals
sigma^2 exactly in the continuum; with a discrete option strip we recover it to
truncation/quadrature error.
"""

import math

import pytest

from quantforge import (
    call_price, put_price, variance_swap_strike, volatility_swap_strike,
)


def _bsm_strip(S0=100, t=1.0, r=0.05, sigma=0.2, step=0.25, hi=600.0):
    Kstar = S0 * math.exp(r * t)
    puts = [k * step for k in range(1, int(Kstar / step))]
    calls = [k * step for k in range(int(Kstar / step), int(hi / step))]
    pp = [put_price(S0, K, t, r, sigma) for K in puts]
    cp = [call_price(S0, K, t, r, sigma) for K in calls]
    return puts, pp, calls, cp, Kstar


@pytest.mark.parametrize("sigma", [0.1, 0.2, 0.4])
def test_fair_variance_matches_flat_vol(sigma):
    S0, t, r = 100, 1.0, 0.05
    puts, pp, calls, cp, Kstar = _bsm_strip(S0, t, r, sigma)
    fair = variance_swap_strike(S0, t, r, puts, pp, calls, cp, split=Kstar)
    # Continuum value is sigma^2; discrete strip recovers it to ~1%.
    assert fair == pytest.approx(sigma * sigma, rel=0.02)


def test_fair_vol_is_sqrt_of_variance():
    S0, t, r, sigma = 100, 1.0, 0.05, 0.3
    puts, pp, calls, cp, Kstar = _bsm_strip(S0, t, r, sigma)
    var = variance_swap_strike(S0, t, r, puts, pp, calls, cp, split=Kstar)
    vol = volatility_swap_strike(S0, t, r, puts, pp, calls, cp, split=Kstar)
    assert vol == pytest.approx(math.sqrt(var), abs=1e-12)
    assert vol == pytest.approx(sigma, rel=0.02)


def test_higher_vol_higher_strike():
    S0, t, r = 100, 1.0, 0.05
    p1, pp1, c1, cp1, ks1 = _bsm_strip(S0, t, r, 0.15)
    p2, pp2, c2, cp2, ks2 = _bsm_strip(S0, t, r, 0.35)
    v1 = variance_swap_strike(S0, t, r, p1, pp1, c1, cp1, split=ks1)
    v2 = variance_swap_strike(S0, t, r, p2, pp2, c2, cp2, split=ks2)
    assert v2 > v1


def test_default_split_is_forward():
    # Passing no split should use the forward and give the same answer.
    S0, t, r, sigma = 100, 1.0, 0.05, 0.2
    puts, pp, calls, cp, Kstar = _bsm_strip(S0, t, r, sigma)
    explicit = variance_swap_strike(S0, t, r, puts, pp, calls, cp, split=Kstar)
    default = variance_swap_strike(S0, t, r, puts, pp, calls, cp)
    assert default == pytest.approx(explicit, abs=1e-9)


def test_rejects_bad_tenor():
    with pytest.raises(ValueError):
        variance_swap_strike(100, 0.0, 0.05, [90, 95], [1, 2], [105, 110], [2, 1])
