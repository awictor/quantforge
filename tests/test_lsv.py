"""LSV leverage-function calibration (Guyon-Henry-Labordere particle method)."""

import math

import pytest

from quantforge import calibrate_lsv_leverage


S0, R = 100.0, 0.0
KAPPA, THETA, XI, RHO, V0 = 1.5, 0.09, 0.5, -0.3, 0.09
EXPIRIES = [0.1, 0.25, 0.5, 1.0]
KGRID = [-0.8 + 0.1 * i for i in range(17)]


def _flat_loc(sigL):
    return lambda K, t: sigL


def test_atm_leverage_is_localvol_over_sqrt_v0():
    # By Gyongy the ATM leverage equals sigma_local / sqrt(E[V|S]); at the very
    # short end E[V|S] = v0, so L_atm = sigL / sqrt(v0).
    sigL = 0.2
    lev, _ = calibrate_lsv_leverage(S0, R, _flat_loc(sigL), KAPPA, THETA, XI,
                                    RHO, V0, EXPIRIES, KGRID, n_paths=8000,
                                    seed=1, sub_steps=10)
    front = lev[EXPIRIES[0]]
    atm = front[min(front, key=abs)]
    assert atm == pytest.approx(sigL / math.sqrt(V0), abs=1e-6)


def test_all_leverage_positive():
    lev, _ = calibrate_lsv_leverage(S0, R, _flat_loc(0.2), KAPPA, THETA, XI,
                                    RHO, V0, EXPIRIES, KGRID, n_paths=6000,
                                    seed=2, sub_steps=8)
    for lt in lev.values():
        assert all(v > 0.0 for v in lt.values())


def test_higher_target_vol_scales_leverage_up():
    lo, _ = calibrate_lsv_leverage(S0, R, _flat_loc(0.15), KAPPA, THETA, XI,
                                   RHO, V0, EXPIRIES, KGRID, n_paths=6000,
                                   seed=3, sub_steps=8)
    hi, _ = calibrate_lsv_leverage(S0, R, _flat_loc(0.30), KAPPA, THETA, XI,
                                   RHO, V0, EXPIRIES, KGRID, n_paths=6000,
                                   seed=3, sub_steps=8)
    lo_atm = lo[EXPIRIES[0]][0.0] if 0.0 in lo[EXPIRIES[0]] else \
        lo[EXPIRIES[0]][min(lo[EXPIRIES[0]], key=abs)]
    hi_atm = hi[EXPIRIES[0]][min(hi[EXPIRIES[0]], key=abs)]
    assert hi_atm > lo_atm


def test_lev_fn_callable_and_positive():
    _, lev_fn = calibrate_lsv_leverage(S0, R, _flat_loc(0.2), KAPPA, THETA, XI,
                                       RHO, V0, EXPIRIES, KGRID, n_paths=6000,
                                       seed=4, sub_steps=8)
    for spot in (85, 100, 120):
        for t in (0.1, 0.5, 1.0):
            v = lev_fn(spot, t)
            assert v > 0.0 and math.isfinite(v)


@pytest.mark.slow
def test_gyongy_identity_holds_on_grid():
    # The calibration sets L = sigma_local / sqrt(condV), so L^2 * condV must
    # equal sigma_local^2 exactly at every calibrated node -- verify via the
    # returned leverage against the flat target.
    sigL = 0.25
    lev, _ = calibrate_lsv_leverage(S0, R, _flat_loc(sigL), KAPPA, THETA, XI,
                                    RHO, V0, EXPIRIES, KGRID, n_paths=10000,
                                    seed=5, sub_steps=15)
    # Front expiry: condV = v0 for all bins, so every L = sigL/sqrt(v0).
    front = lev[EXPIRIES[0]]
    for L in front.values():
        assert L == pytest.approx(sigL / math.sqrt(V0), abs=1e-6)
