"""Pearson correlation significance test and Fisher-z confidence interval."""

import math
import random

import pytest

from quantforge import pearson_r, pearson_correlation_test
from quantforge.mathfns import norm_ppf


def test_pearson_matches_manual():
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 5, 4, 5]
    mx, my = 3.0, 4.0
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    sxx = sum((xi - mx) ** 2 for xi in x)
    syy = sum((yi - my) ** 2 for yi in y)
    assert abs(pearson_r(x, y) - sxy / math.sqrt(sxx * syy)) < 1e-12


def test_t_statistic_identity():
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 5, 4, 5]
    res = pearson_correlation_test(x, y)
    r, df = res["r"], res["df"]
    assert abs(res["t_stat"] - r * math.sqrt(df / (1 - r * r))) < 1e-10


def test_perfect_correlation():
    res = pearson_correlation_test([1, 2, 3, 4], [2, 4, 6, 8])
    assert abs(res["r"] - 1.0) < 1e-12
    assert res["p_value"] == 0.0
    assert res["conf_int"] == [res["r"], res["r"]]


def test_perfect_negative_correlation():
    res = pearson_correlation_test([1, 2, 3, 4], [8, 6, 4, 2])
    assert abs(res["r"] + 1.0) < 1e-12


def test_confidence_interval_contains_r_for_strong_signal():
    x = list(range(20))
    y = [2 * xi + 1 for xi in x]
    res = pearson_correlation_test(x, y)
    lo, hi = res["conf_int"]
    assert lo <= res["r"] <= hi + 1e-12


def test_fisher_interval_coverage():
    rho, n = 0.6, 30
    rng = random.Random(7)
    hits, trials = 0, 3000
    for _ in range(trials):
        xs, ys = [], []
        for _ in range(n):
            z1, z2 = rng.gauss(0, 1), rng.gauss(0, 1)
            xs.append(z1)
            ys.append(rho * z1 + math.sqrt(1 - rho * rho) * z2)
        lo, hi = pearson_correlation_test(xs, ys)["conf_int"]
        if lo <= rho <= hi:
            hits += 1
    assert abs(hits / trials - 0.95) < 0.02


def test_higher_confidence_widens_interval():
    x = list(range(15))
    y = [xi + (xi % 3) for xi in x]
    w90 = pearson_correlation_test(x, y, 0.90)["conf_int"]
    w99 = pearson_correlation_test(x, y, 0.99)["conf_int"]
    assert (w99[1] - w99[0]) > (w90[1] - w90[0])


def test_strong_correlation_significant_noise_not():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(80)]
    y_signal = [2 * xi + rng.gauss(0, 0.1) for xi in x]
    y_noise = [rng.gauss(0, 1) for _ in range(80)]
    assert pearson_correlation_test(x, y_signal)["p_value"] < 1e-20
    assert pearson_correlation_test(x, y_noise)["p_value"] > 0.05


def test_validation():
    with pytest.raises(ValueError):
        pearson_r([1.0], [1.0])
    with pytest.raises(ValueError):
        pearson_r([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])   # zero variance in x
    with pytest.raises(ValueError):
        pearson_correlation_test([1.0, 2.0], [1.0, 2.0])   # n < 3
    with pytest.raises(ValueError):
        pearson_correlation_test([1, 2, 3], [1, 2, 3], confidence=1.0)
