"""Double-Heston calibration to an implied-vol surface."""

import pytest

from quantforge import (
    OptionType,
    calibrate_double_heston,
    double_heston_price,
)
from quantforge.implied import implied_volatility


S, R = 100.0, 0.02
TRUE = (0.04, 3.0, 0.04, 0.5, -0.7, 0.03, 0.4, 0.05, 0.25, -0.4)
EXPIRIES = [0.25, 0.5, 1.0, 2.0]
STRIKES = [85, 92, 100, 108, 116]


def _surface(params=TRUE):
    quotes = []
    for t in EXPIRIES:
        for K in STRIKES:
            c = double_heston_price(S, K, t, R, *params, OptionType.CALL)
            iv = implied_volatility(c, S, K, t, R, OptionType.CALL, b=R)
            quotes.append((t, K, iv))
    return quotes


@pytest.mark.slow
def test_calibration_fits_surface():
    quotes = _surface()
    params, rmse = calibrate_double_heston(S, R, quotes, max_iter=4000)
    assert rmse < 1e-3
    # Spot-check that the fitted params reprice a few quotes.
    for t, K, mv in quotes[:5]:
        c = double_heston_price(S, K, t, R, *params, OptionType.CALL)
        iv = implied_volatility(c, S, K, t, R, OptionType.CALL, b=R)
        assert iv == pytest.approx(mv, abs=2e-3)


@pytest.mark.slow
def test_calibrated_params_are_valid():
    quotes = _surface()
    params, _ = calibrate_double_heston(S, R, quotes, max_iter=3000)
    v01, k1, th1, xi1, rho1, v02, k2, th2, xi2, rho2 = params
    for v in (v01, th1, xi1, v02, th2, xi2):
        assert v > 0.0
    for k in (k1, k2):
        assert k > 0.0
    for rho in (rho1, rho2):
        assert -1.0 < rho < 1.0


def test_requires_enough_quotes():
    with pytest.raises(ValueError):
        calibrate_double_heston(S, R, [(1.0, 100, 0.2)])
