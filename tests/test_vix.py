"""CBOE VIX-style fair volatility index from a chain or smile."""

import math

import pytest

from quantforge import vix_from_chain, vix_from_smile
from quantforge.bsm import call_price, put_price


S0, T, R = 100.0, 30.0 / 365.0, 0.03


@pytest.mark.parametrize("sig", [0.15, 0.2, 0.35])
def test_flat_smile_gives_100_sigma(sig):
    _var, vix = vix_from_smile(S0, T, R, lambda K: sig,
                               n_strikes=301, width=8.0)
    assert vix == pytest.approx(100.0 * sig, abs=0.1)


def test_downward_skew_lifts_vix_above_atm():
    def smile(K):
        return max(0.05, 0.2 + 0.15 * math.log(S0 / K))
    _var, vix = vix_from_smile(S0, T, R, smile, n_strikes=301, width=8.0)
    assert vix > 20.0


def test_chain_matches_flat_vol():
    sig = 0.2
    F = S0 * math.exp(R * T)
    strikes = [F * math.exp((-8 + 16 * i / 200) * sig * math.sqrt(T))
               for i in range(201)]
    prices = [put_price(S0, K, T, R, sig, b=R) if K < F
              else call_price(S0, K, T, R, sig, b=R) for K in strikes]
    _var, vix = vix_from_chain(strikes, prices, F, T, R)
    assert vix == pytest.approx(20.0, abs=0.1)


def test_variance_is_vix_squared():
    var, vix = vix_from_smile(S0, T, R, lambda K: 0.25, n_strikes=201)
    assert vix == pytest.approx(100.0 * math.sqrt(var))


def test_short_chain_raises():
    with pytest.raises(ValueError):
        vix_from_chain([100, 110], [5, 3], 100, T, R)


def test_bad_tenor_raises():
    with pytest.raises(ValueError):
        vix_from_smile(S0, 0.0, R, lambda K: 0.2)
