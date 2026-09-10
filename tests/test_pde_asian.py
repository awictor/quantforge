"""Arithmetic-Asian augmented-state PDE vs closed form and Monte Carlo."""

import math

import pytest

from quantforge import (
    OptionType,
    asian_pde_price,
    arithmetic_asian,
)


S, K, T, R, SIGMA = 100.0, 100.0, 1.0, 0.05, 0.2


def test_zero_strike_is_discounted_average():
    # A K->0 average-price call is the discounted expected average.
    b = R
    exact = math.exp(-R * T) * S * (math.exp(b * T) - 1.0) / (b * T)
    pde = asian_pde_price(S, 1e-6, T, R, SIGMA, OptionType.CALL,
                          n_s=150, n_i=300, n_time=150)
    assert pde == pytest.approx(exact, abs=0.05)


@pytest.mark.slow
@pytest.mark.parametrize("S0,Kk,sig", [(100, 100, 0.2), (100, 90, 0.3),
                                       (100, 110, 0.25)])
def test_close_to_turnbull_wakeman(S0, Kk, sig):
    tw = arithmetic_asian(S0, Kk, T, R, sig, OptionType.CALL)
    pde = asian_pde_price(S0, Kk, T, R, sig, OptionType.CALL,
                          n_s=150, n_i=300, n_time=150)
    # Both are approximations of the same continuous average; agree to ~0.25.
    assert pde == pytest.approx(tw, abs=0.25)


@pytest.mark.slow
def test_converges_toward_monte_carlo():
    from quantforge.montecarlo import arithmetic_asian_mc
    mc = arithmetic_asian_mc(S, K, T, R, SIGMA, OptionType.CALL,
                             n_steps=200, n_paths=200_000, seed=1)
    e_coarse = abs(asian_pde_price(S, K, T, R, SIGMA, n_s=120, n_i=200,
                                   n_time=120) - mc.price)
    e_fine = abs(asian_pde_price(S, K, T, R, SIGMA, n_s=200, n_i=500,
                                 n_time=200) - mc.price)
    assert e_fine < e_coarse


def test_put_positive_and_below_strike_disc():
    p = asian_pde_price(S, 110, T, R, SIGMA, OptionType.PUT,
                        n_s=120, n_i=250, n_time=120)
    assert 0.0 < p < 110 * math.exp(-R * T)


def test_zero_time_is_intrinsic():
    assert asian_pde_price(120, 100, 0.0, R, SIGMA,
                           OptionType.CALL) == pytest.approx(20.0)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        asian_pde_price(-1, 100, T, R, SIGMA)
