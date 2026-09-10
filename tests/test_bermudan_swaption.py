"""Bermudan swaption via Longstaff-Schwartz on the G2++ state."""

import math

import pytest

from quantforge import bermudan_swaption_g2pp
from quantforge.bermudan_swaption import _swap_value, _simulate_state


R0 = 0.03
A, B, SIGMA, ETA, RHO = 0.1, 0.3, 0.01, 0.008, -0.3


def _P0(t):
    return math.exp(-R0 * t)


def test_price_positive():
    v = bermudan_swaption_g2pp(_P0, [1.0, 2.0, 3.0], 0.03, A, B, SIGMA, ETA, RHO,
                               payer=True, n_paths=20000, seed=1)
    assert v > 0.0 and math.isfinite(v)


@pytest.mark.slow
def test_bermudan_at_least_single_exercise():
    # More exercise opportunities cannot be worth less than one.
    eur = bermudan_swaption_g2pp(_P0, [1.0], 0.03, A, B, SIGMA, ETA, RHO,
                                 payer=True, n_paths=40000, seed=2)
    berm = bermudan_swaption_g2pp(_P0, [1.0, 2.0, 3.0, 4.0], 0.03,
                                  A, B, SIGMA, ETA, RHO,
                                  payer=True, n_paths=40000, seed=2)
    assert berm >= eur - 1e-3


def test_payer_and_receiver_both_positive():
    payer = bermudan_swaption_g2pp(_P0, [1.0, 2.0, 3.0], 0.03,
                                   A, B, SIGMA, ETA, RHO, payer=True,
                                   n_paths=6000, seed=3)
    recv = bermudan_swaption_g2pp(_P0, [1.0, 2.0, 3.0], 0.03,
                                  A, B, SIGMA, ETA, RHO, payer=False,
                                  n_paths=6000, seed=3)
    assert payer > 0.0 and recv > 0.0


def test_higher_strike_lowers_payer_value():
    lo = bermudan_swaption_g2pp(_P0, [1.0, 2.0, 3.0], 0.02, A, B, SIGMA, ETA, RHO,
                                payer=True, n_paths=6000, seed=4)
    hi = bermudan_swaption_g2pp(_P0, [1.0, 2.0, 3.0], 0.06, A, B, SIGMA, ETA, RHO,
                                payer=True, n_paths=6000, seed=4)
    assert lo > hi  # a payer is worth less as the fixed rate rises


def test_swap_value_payer_receiver_antisymmetric():
    sv_p = _swap_value(_P0, 0.001, -0.001, A, B, SIGMA, ETA, RHO, 1.0, 0.03,
                       [2, 3, 4], True)
    sv_r = _swap_value(_P0, 0.001, -0.001, A, B, SIGMA, ETA, RHO, 1.0, 0.03,
                       [2, 3, 4], False)
    assert sv_p == pytest.approx(-sv_r, abs=1e-12)


def test_state_simulation_shape_and_mean():
    import random
    rng = random.Random(0)
    paths = _simulate_state(A, B, SIGMA, ETA, RHO, [1.0, 2.0], 5000, rng)
    assert len(paths) == 5000 and len(paths[0]) == 2
    # OU factors have zero mean under the (approximate) simulated measure.
    mean_x = sum(p[0][0] for p in paths) / len(paths)
    assert abs(mean_x) < 0.01


def test_requires_exercise_dates():
    with pytest.raises(ValueError):
        bermudan_swaption_g2pp(_P0, [], 0.03, A, B, SIGMA, ETA, RHO)
