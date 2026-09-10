"""SSVI local-vol Monte Carlo closes the calibrate -> local-vol -> reprice loop."""

import math

import pytest

from quantforge import (
    SSVIParams,
    call_price,
    ssvi_local_vol_fn,
    ssvi_reprice_mc,
)


# A surface whose expiries reach down near the origin, so a Monte Carlo path
# integrating tau from 0 sees the true short-maturity local vol.
PARAMS = SSVIParams(
    rho=-0.4, eta=1.0, gamma=0.5,
    thetas={0.05: 0.045 * 0.05, 0.25: 0.045 * 0.25,
            0.5: 0.045 * 0.5, 1.0: 0.045},
)
S0, R = 100.0, 0.0


@pytest.mark.slow
@pytest.mark.parametrize("K", [80, 90, 100, 110, 120])
def test_local_vol_mc_reprices_ssvi_smile(K):
    t = 1.0
    F = S0 * math.exp(R * t)
    k = math.log(K / F)
    bs = call_price(S0, K, t, R, PARAMS.implied_vol(k, t), b=R)
    mc = ssvi_reprice_mc(PARAMS, S0, K, t, R, n_steps=120, n_paths=100_000, seed=1)
    # The local-vol surface must reproduce the SSVI implied price to within MC
    # error (a small Euler bias aside).
    assert abs(mc.price - bs) < 3.0 * mc.std_error + 0.05


def test_local_vol_fn_atm_matches_implied_short_time():
    # As tau -> 0 the local vol at the money approaches the ATM implied vol.
    lv = ssvi_local_vol_fn(PARAMS, S0, R)
    atm_implied = PARAMS.implied_vol(0.0, 0.05)
    assert lv(S0, 0.05) == pytest.approx(atm_implied, abs=2e-2)


def test_local_vol_fn_is_callable_and_positive():
    lv = ssvi_local_vol_fn(PARAMS, S0, R)
    for spot in (80.0, 100.0, 125.0):
        for tau in (0.1, 0.5, 1.0):
            v = lv(spot, tau)
            assert v > 0.0 and math.isfinite(v)


def test_from_params_rejects_beyond_longest_expiry():
    from quantforge import ssvi_local_vol_from_params
    with pytest.raises(ValueError):
        ssvi_local_vol_from_params(PARAMS, 0.0, 5.0)
