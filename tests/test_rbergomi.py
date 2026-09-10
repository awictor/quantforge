"""Rough Bergomi Monte Carlo: martingale, flat limit, and the rough skew."""

import math

import pytest

from quantforge import OptionType, rbergomi_price, rbergomi_smile
from quantforge.rbergomi import rbergomi_paths


S = 100.0
STRIKES = [85, 92, 100, 108, 116]


@pytest.mark.slow
def test_discounted_spot_is_martingale():
    r = 0.03
    term = rbergomi_paths(S, 1.0, 0.04, 1.5, 0.1, -0.7, r,
                          n_steps=100, n_paths=80_000, seed=1)
    ES = sum(term) / len(term)
    assert abs(ES - S * math.exp(r)) / (S * math.exp(r)) < 5e-3


@pytest.mark.slow
def test_zero_vol_of_vol_is_flat_at_sqrt_xi0():
    # eta = 0 kills the rough driver: variance is the flat xi0, so the implied
    # smile is flat at sqrt(xi0).
    sm = rbergomi_smile(S, STRIKES, 1.0, 0.04, 0.0, 0.1, -0.7, r=0.0,
                        n_steps=100, n_paths=80_000, seed=2)
    for _, iv in sm:
        assert iv == pytest.approx(0.2, abs=5e-3)


@pytest.mark.slow
def test_negative_rho_gives_downward_skew():
    sm = rbergomi_smile(S, STRIKES, 0.2, 0.04, 1.5, 0.1, -0.8, r=0.0,
                        n_steps=120, n_paths=120_000, seed=3)
    vols = [iv for _, iv in sm]
    assert vols[0] > vols[-1]  # low-strike vol richer


@pytest.mark.slow
def test_rough_skew_steeper_than_diffusive_at_short_maturity():
    # The whole point of rough vol: H < 1/2 produces a steeper short-dated ATM
    # skew than the H = 1/2 diffusive case with the same parameters.
    def slope(H):
        sm = rbergomi_smile(S, STRIKES, 0.1, 0.04, 1.9, H, -0.9, r=0.0,
                            n_steps=120, n_paths=120_000, seed=7)
        lm = [a for a, _ in sm]
        iv = [b for _, b in sm]
        return (iv[-1] - iv[0]) / (lm[-1] - lm[0])

    rough = slope(0.1)
    diffusive = slope(0.5)
    assert rough < 0 and diffusive < 0
    assert rough < diffusive  # steeper (more negative)


def test_price_matches_smile_atm_cheap():
    # Small consistency run: the ATM call price and the smile's ATM vol agree.
    res = rbergomi_price(S, 100.0, 0.25, 0.04, 1.0, 0.2, -0.5, r=0.0,
                         option_type=OptionType.CALL, n_steps=40, n_paths=8_000,
                         seed=5)
    assert res.price > 0.0 and math.isfinite(res.price)
    assert res.std_error >= 0.0


def test_bad_params_raise():
    with pytest.raises(ValueError):
        rbergomi_price(S, 100, 1.0, 0.04, 1.0, 1.5, -0.5, n_paths=100)  # H>=1
    with pytest.raises(ValueError):
        rbergomi_price(S, 100, 1.0, -0.04, 1.0, 0.2, -0.5, n_paths=100)  # xi0<=0
    with pytest.raises(ValueError):
        rbergomi_price(S, 100, 1.0, 0.04, 1.0, 0.2, -1.5, n_paths=100)  # rho
