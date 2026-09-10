"""Heston 5-parameter calibration to an implied-vol surface."""

import pytest

from quantforge import OptionType, calibrate_heston, heston_price
from quantforge.implied import implied_volatility


S, R = 100.0, 0.02
TRUE = (0.04, 1.5, 0.05, 0.4, -0.6)
EXPIRIES = [0.25, 0.5, 1.0, 2.0]
STRIKES = [85, 92, 100, 108, 116]


def _surface(params=TRUE):
    quotes = []
    for t in EXPIRIES:
        for K in STRIKES:
            c = heston_price(S, K, t, R, *params, OptionType.CALL)
            iv = implied_volatility(c, S, K, t, R, OptionType.CALL, b=R)
            quotes.append((t, K, iv))
    return quotes


@pytest.mark.slow
def test_recovers_synthetic_parameters():
    params, rmse = calibrate_heston(S, R, _surface(), max_iter=4000)
    assert rmse < 1e-3
    for got, want in zip(params, TRUE):
        assert got == pytest.approx(want, abs=0.05)


@pytest.mark.slow
def test_feller_penalty_pushes_toward_condition():
    # A large Feller weight should keep 2 kappa theta from falling far below xi^2.
    params, _ = calibrate_heston(S, R, _surface(), feller_weight=10.0,
                                 max_iter=4000)
    v0, kappa, theta, xi, rho = params
    assert 2.0 * kappa * theta >= xi * xi - 1e-3


def test_calibrated_params_valid():
    # Validity of the reparametrization, not fit accuracy -- a short run suffices.
    params, _ = calibrate_heston(S, R, _surface(), max_iter=300)
    v0, kappa, theta, xi, rho = params
    assert v0 > 0 and kappa > 0 and theta > 0 and xi > 0
    assert -1.0 < rho < 1.0


def test_requires_enough_quotes():
    with pytest.raises(ValueError):
        calibrate_heston(S, R, [(1.0, 100, 0.2), (1.0, 110, 0.19)])
