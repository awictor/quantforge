"""GARCH term (average) volatility over a horizon (volatility.garch_term_variance)."""

import math

import pytest

from quantforge import GarchParams, garch_forecast, garch_term_variance


P = GarchParams(omega=2e-6, alpha=0.08, beta=0.90)
LR = P.long_run_variance
H = LR * 1.5   # elevated starting variance
R = 0.01


def _brute(n):
    total = 0.0
    for k in range(1, n + 1):
        total += garch_forecast(P, R, H, k, 252) ** 2 / 252  # de-annualize
    return math.sqrt(total / n * 252)


@pytest.mark.parametrize("n", [1, 5, 20, 60])
def test_matches_average_of_point_forecasts(n):
    assert garch_term_variance(P, R, H, n, 252) == pytest.approx(_brute(n), abs=1e-9)


def test_horizon_one_equals_point_forecast():
    assert garch_term_variance(P, R, H, 1, 252) == pytest.approx(
        garch_forecast(P, R, H, 1, 252), abs=1e-12)


def test_converges_to_long_run_vol():
    lr_vol = math.sqrt(LR * 252)
    assert garch_term_variance(P, R, H, 100000, 252) == pytest.approx(lr_vol, abs=1e-4)


def test_elevated_variance_term_decreases_with_horizon():
    short = garch_term_variance(P, R, H, 5, 252)
    long = garch_term_variance(P, R, H, 120, 252)
    assert short > long > math.sqrt(LR * 252)


def test_unit_root_uses_one_step():
    p = GarchParams(omega=0.0, alpha=0.1, beta=0.9)  # persistence == 1
    h1 = garch_forecast(p, R, H, 1, 252)
    assert garch_term_variance(p, R, H, 10, 252) == pytest.approx(h1, abs=1e-12)


def test_validation():
    with pytest.raises(ValueError):
        garch_term_variance(P, R, H, 0, 252)
