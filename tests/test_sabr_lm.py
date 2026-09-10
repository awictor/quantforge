"""Levenberg-Marquardt SABR fit (analytic Jacobian) vs Nelder-Mead."""

import pytest

from quantforge import (
    SABRParams,
    sabr_vol,
    calibrate_sabr,
    calibrate_sabr_lm,
)


F, T, BETA = 100.0, 1.0, 0.5
STRIKES = [70, 85, 100, 115, 130, 150]
TRUE = SABRParams(alpha=2.0, beta=BETA, rho=-0.35, nu=0.45)


def _smile(params=TRUE):
    return [sabr_vol(F, K, T, params.alpha, BETA, params.rho, params.nu)
            for K in STRIKES]


def test_recovers_synthetic_params_fast():
    vols = _smile()
    params, rmse, n_iter = calibrate_sabr_lm(F, T, STRIKES, vols, beta=BETA)
    assert params.alpha == pytest.approx(TRUE.alpha, abs=1e-4)
    assert params.rho == pytest.approx(TRUE.rho, abs=1e-4)
    assert params.nu == pytest.approx(TRUE.nu, abs=1e-4)
    assert rmse < 1e-6
    assert n_iter <= 20  # converges in a handful of Gauss-Newton steps


def test_matches_nelder_mead_on_noisy_market():
    import random
    rng = random.Random(1)
    vols = [v * (1.0 + rng.uniform(-0.01, 0.01)) for v in _smile()]
    lm_params, lm_rmse, _ = calibrate_sabr_lm(F, T, STRIKES, vols, beta=BETA)
    nm_params, nm_rmse = calibrate_sabr(F, T, STRIKES, vols, beta=BETA)
    # Both minimisers should land on essentially the same optimum.
    assert lm_rmse == pytest.approx(nm_rmse, rel=0.02)
    assert lm_params.alpha == pytest.approx(nm_params.alpha, rel=0.02)
    assert lm_params.nu == pytest.approx(nm_params.nu, rel=0.05)


def test_result_refits_the_input_smile():
    vols = _smile()
    params, _, _ = calibrate_sabr_lm(F, T, STRIKES, vols, beta=BETA)
    for K, v in zip(STRIKES, vols):
        model = sabr_vol(F, K, T, params.alpha, BETA, params.rho, params.nu)
        assert model == pytest.approx(v, abs=1e-5)


def test_respects_parameter_bounds():
    vols = _smile()
    params, _, _ = calibrate_sabr_lm(F, T, STRIKES, vols, beta=BETA)
    assert params.alpha > 0.0
    assert params.nu >= 0.0
    assert -1.0 < params.rho < 1.0


def test_requires_three_points():
    with pytest.raises(ValueError):
        calibrate_sabr_lm(F, T, [90, 110], [0.2, 0.19], beta=BETA)


def test_beta_one_lognormal_smile():
    true = SABRParams(alpha=0.2, beta=1.0, rho=-0.4, nu=0.5)
    vols = [sabr_vol(F, K, T, true.alpha, 1.0, true.rho, true.nu)
            for K in STRIKES]
    params, rmse, _ = calibrate_sabr_lm(F, T, STRIKES, vols, beta=1.0)
    assert rmse < 1e-6
    assert params.rho == pytest.approx(true.rho, abs=1e-3)
