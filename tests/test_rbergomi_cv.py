"""Conditional (turbocharged) rough Bergomi estimator: agreement + variance cut."""

import math

import pytest

from quantforge import (
    OptionType,
    rbergomi_price,
    rbergomi_price_cv,
    rbergomi_smile,
    rbergomi_smile_cv,
)


S = 100.0
STRIKES = [85, 92, 100, 108, 116]
PARAMS = dict(xi0=0.04, eta=1.5, H=0.1, rho=-0.7, r=0.0)


@pytest.mark.slow
def test_cv_agrees_with_plain_and_cuts_se():
    t = 0.5
    plain = rbergomi_price(S, 100.0, t, **PARAMS, option_type=OptionType.CALL,
                           n_steps=80, n_paths=60_000, seed=1)
    cv = rbergomi_price_cv(S, 100.0, t, **PARAMS,
                           n_steps=80, n_paths=60_000, seed=1)
    # Same price within combined error, but a materially smaller standard error.
    assert abs(plain.price - cv.price) < 3.0 * (plain.std_error + cv.std_error)
    assert cv.std_error < plain.std_error


@pytest.mark.slow
def test_se_reduction_grows_as_rho_shrinks():
    # The orthogonal-noise fraction integrated out is 1 - rho^2, so the variance
    # reduction should improve as |rho| decreases.
    t = 0.5

    def ratio(rho):
        p = dict(PARAMS)
        p["rho"] = rho
        a = rbergomi_price(S, 100.0, t, **p, option_type=OptionType.CALL,
                           n_steps=80, n_paths=40_000, seed=2)
        b = rbergomi_price_cv(S, 100.0, t, **p,
                              n_steps=80, n_paths=40_000, seed=2)
        return a.std_error / b.std_error

    assert ratio(-0.1) > ratio(-0.9)


@pytest.mark.slow
def test_cv_smile_matches_plain_smile():
    t = 0.5
    cv = rbergomi_smile_cv(S, STRIKES, t, **PARAMS,
                           n_steps=80, n_paths=60_000, seed=5)
    pl = rbergomi_smile(S, STRIKES, t, **PARAMS,
                        n_steps=80, n_paths=60_000, seed=5)
    for (_, c), (_, p) in zip(cv, pl):
        assert c == pytest.approx(p, abs=5e-3)


def test_fully_correlated_branch():
    # rho = -1 leaves no residual noise; the conditional price is well-defined.
    p = dict(PARAMS)
    p["rho"] = -1.0
    res = rbergomi_price_cv(S, 100.0, 0.5, **p,
                            n_steps=50, n_paths=4_000, seed=3)
    assert res.price > 0.0 and math.isfinite(res.price)


def test_bad_params_raise():
    with pytest.raises(ValueError):
        rbergomi_price_cv(S, 100.0, 0.5, xi0=0.04, eta=1.0, H=1.5, rho=-0.5,
                          n_paths=100)
    with pytest.raises(ValueError):
        rbergomi_price_cv(S, -1.0, 0.5, xi0=0.04, eta=1.0, H=0.2, rho=-0.5,
                          n_paths=100)
