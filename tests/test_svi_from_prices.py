"""Raw SVI calibration directly from market call prices."""

import math
import random

import pytest

from quantforge import SVIParams, calibrate_svi_from_prices, call_price


F, T, R = 100.0, 1.0, 0.03
S = F * math.exp(-R * T)
TRUE = SVIParams(a=0.04, b=0.4, rho=-0.4, m=0.0, s=0.1)
STRIKES = [70, 80, 90, 100, 110, 120, 130]


def _prices(params=TRUE):
    out = []
    for K in STRIKES:
        k = math.log(K / F)
        iv = math.sqrt(params.total_variance(k) / T)
        out.append(call_price(S, K, T, R, iv, b=R))
    return out


def test_recovers_svi_from_clean_prices():
    params, iv_rmse, px_rmse = calibrate_svi_from_prices(F, T, R, STRIKES,
                                                         _prices())
    assert iv_rmse < 1e-4
    assert px_rmse < 1e-3
    assert params.a == pytest.approx(TRUE.a, abs=1e-3)
    assert params.rho == pytest.approx(TRUE.rho, abs=1e-2)


def test_fits_noisy_prices():
    rng = random.Random(1)
    noisy = [c * (1 + rng.uniform(-0.005, 0.005)) for c in _prices()]
    _params, _iv, px_rmse = calibrate_svi_from_prices(F, T, R, STRIKES, noisy)
    assert px_rmse < 0.1


def test_vega_and_unweighted_both_fit():
    prices = _prices()
    p1, _, r1 = calibrate_svi_from_prices(F, T, R, STRIKES, prices,
                                          vega_weighted=True)
    p0, _, r0 = calibrate_svi_from_prices(F, T, R, STRIKES, prices,
                                          vega_weighted=False)
    assert r1 < 1e-3 and r0 < 1e-3


def test_reprices_within_price_rmse():
    params, _, px_rmse = calibrate_svi_from_prices(F, T, R, STRIKES, _prices())
    for K, c in zip(STRIKES, _prices()):
        k = math.log(K / F)
        iv = math.sqrt(params.total_variance(k) / T)
        model = call_price(S, K, T, R, iv, b=R)
        assert abs(model - c) < 5e-3


def test_requires_enough_quotes():
    with pytest.raises(ValueError):
        calibrate_svi_from_prices(F, T, R, [90, 100, 110], [12, 6, 3])
