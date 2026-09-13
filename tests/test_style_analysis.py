"""Sharpe returns-based style analysis."""

import random

import pytest

from quantforge import style_analysis


def _indices(seed=1, n=500):
    rng = random.Random(seed)
    a = [rng.gauss(0.0005, 0.012) for _ in range(n)]
    b = [rng.gauss(0.0003, 0.008) for _ in range(n)]
    c = [rng.gauss(0.0002, 0.02) for _ in range(n)]
    return a, b, c


def test_recovers_known_blend():
    a, b, c = _indices()
    n = len(a)
    rng = random.Random(9)
    fund = [0.6 * a[i] + 0.4 * b[i] + rng.gauss(0, 0.0005) for i in range(n)]
    res = style_analysis(fund, [a, b, c])
    w = res["weights"]
    assert abs(w[0] - 0.6) < 0.05
    assert abs(w[1] - 0.4) < 0.05
    assert w[2] < 0.05
    assert res["r_squared"] > 0.95


def test_weights_on_simplex():
    a, b, c = _indices()
    n = len(a)
    fund = [0.5 * a[i] + 0.5 * c[i] for i in range(n)]
    res = style_analysis(fund, [a, b, c])
    assert abs(sum(res["weights"]) - 1.0) < 1e-6
    assert all(w >= -1e-12 for w in res["weights"])


def test_exact_blend_full_r_squared():
    a, b, c = _indices()
    n = len(a)
    fund = [0.7 * a[i] + 0.3 * c[i] for i in range(n)]
    res = style_analysis(fund, [a, b, c])
    assert res["r_squared"] > 0.999
    assert res["tracking_error"] < 1e-4


def test_single_index():
    a, _, _ = _indices()
    res = style_analysis(a, [a])
    assert abs(res["weights"][0] - 1.0) < 1e-9
    assert res["r_squared"] > 0.999


def test_validation():
    with pytest.raises(ValueError):
        style_analysis([0.01, 0.02], [])
    with pytest.raises(ValueError):
        style_analysis([0.01, 0.02], [[0.01]])       # length mismatch
