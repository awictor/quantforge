"""AR(p) fitting by Yule-Walker and forecasting."""

import random

import pytest

from quantforge import fit_ar_yule_walker, ar_forecast


def test_recovers_ar1_parameters():
    random.seed(1)
    phi, c = 0.6, 2.0
    y = [5.0]
    for _ in range(20000):
        y.append(c + phi * y[-1] + random.gauss(0, 1))
    m = fit_ar_yule_walker(y[1:], 1)
    assert abs(m["coefficients"][0] - phi) < 0.03
    assert abs(m["intercept"] - c) < 0.2
    assert abs(m["noise_variance"] - 1.0) < 0.1


def test_recovers_ar2_parameters():
    random.seed(2)
    a1, a2 = 0.5, -0.3
    y = [0.0, 0.0]
    for _ in range(20000):
        y.append(a1 * y[-1] + a2 * y[-2] + random.gauss(0, 1))
    m = fit_ar_yule_walker(y[2:], 2)
    assert abs(m["coefficients"][0] - a1) < 0.03
    assert abs(m["coefficients"][1] - a2) < 0.03


def test_forecast_mean_reverts():
    random.seed(1)
    y = [5.0]
    for _ in range(20000):
        y.append(2.0 + 0.6 * y[-1] + random.gauss(0, 1))
    m = fit_ar_yule_walker(y[1:], 1)
    f = ar_forecast(m, [8.0], steps=5)
    assert f[0] < 8.0
    assert abs(f[-1] - m["mean"]) < abs(f[0] - m["mean"])


def test_one_step_matches_recursion():
    random.seed(1)
    y = [5.0]
    for _ in range(5000):
        y.append(2.0 + 0.6 * y[-1] + random.gauss(0, 1))
    m = fit_ar_yule_walker(y[1:], 1)
    expected = m["intercept"] + m["coefficients"][0] * 10.0
    assert abs(ar_forecast(m, [10.0], 1)[0] - expected) < 1e-9


def test_validation():
    rng = random.Random(1)
    y = [rng.gauss(0, 1) for _ in range(100)]
    with pytest.raises(ValueError):
        fit_ar_yule_walker(y, 0)
    with pytest.raises(ValueError):
        fit_ar_yule_walker([1.0, 2.0], 3)     # too short
    m = fit_ar_yule_walker(y, 1)
    with pytest.raises(ValueError):
        ar_forecast(m, [], 1)                 # history shorter than order
