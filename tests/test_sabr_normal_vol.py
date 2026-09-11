"""Hagan normal (Bachelier) SABR implied volatility (sabr_normal_vol)."""

import pytest

from quantforge import sabr_normal_vol, sabr_vol
from quantforge.bsm import price as bsm_price
from quantforge.bachelier import bachelier_implied_vol


F, T = 100.0, 1.0
ALPHA, BETA, RHO, NU = 0.2, 0.5, -0.3, 0.4
STRIKES = [80.0, 90.0, 100.0, 110.0, 120.0]


@pytest.mark.parametrize("K", STRIKES)
def test_matches_black_to_bachelier_roundtrip(K):
    # Reference normal vol: SABR-Black vol -> BS price -> Bachelier implied vol.
    black = sabr_vol(F, K, T, ALPHA, BETA, RHO, NU)
    pr = bsm_price(F, K, T, 0.0, black, "call", b=0.0)
    ref = bachelier_implied_vol(pr, F, K, T, 0.0, "call")
    assert sabr_normal_vol(F, K, T, ALPHA, BETA, RHO, NU) == pytest.approx(ref, abs=5e-3)


def test_atm_limit_continuous():
    atm = sabr_normal_vol(F, F, T, ALPHA, BETA, RHO, NU)
    near = sabr_normal_vol(F, F * (1 + 1e-7), T, ALPHA, BETA, RHO, NU)
    assert atm == pytest.approx(near, rel=1e-4)


def test_beta_zero_no_vol_of_vol_is_alpha():
    # beta = 0, nu = 0: ATM normal vol is exactly alpha.
    v = sabr_normal_vol(F, F, T, 0.15, 0.0, 0.0, 0.0)
    assert v == pytest.approx(0.15, rel=1e-9)


def test_positive_and_smile_shape():
    vols = [sabr_normal_vol(F, K, T, ALPHA, BETA, RHO, NU) for K in STRIKES]
    assert all(v > 0 for v in vols)
    atm_idx = STRIKES.index(100.0)
    assert vols[0] > vols[atm_idx]
    assert vols[-1] > vols[atm_idx]


def test_validation():
    with pytest.raises(ValueError):
        sabr_normal_vol(-1.0, 100.0, T, ALPHA, BETA, RHO, NU)
    with pytest.raises(ValueError):
        sabr_normal_vol(F, 100.0, T, -0.1, BETA, RHO, NU)
    with pytest.raises(ValueError):
        sabr_normal_vol(F, 100.0, 0.0, ALPHA, BETA, RHO, NU)
