"""Forecast combination weights."""

import random

import pytest

from quantforge import (simple_average_forecast, inverse_mse_weights,
                        optimal_combination_weights, combine_forecasts)


def _mse(e):
    return sum(x * x for x in e) / len(e)


def test_simple_average():
    f = simple_average_forecast([[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]])
    assert f == [2.0, 2.0, 2.0]


def test_weights_sum_to_one():
    rng = random.Random(1)
    e1 = [rng.gauss(0, 1) for _ in range(500)]
    e2 = [rng.gauss(0, 2) for _ in range(500)]
    assert abs(sum(inverse_mse_weights([e1, e2])) - 1.0) < 1e-12
    assert abs(sum(optimal_combination_weights([e1, e2])) - 1.0) < 1e-12


def test_inverse_mse_favors_accurate_model():
    rng = random.Random(2)
    accurate = [rng.gauss(0, 0.5) for _ in range(500)]
    noisy = [rng.gauss(0, 2.0) for _ in range(500)]
    w = inverse_mse_weights([accurate, noisy])
    assert w[0] > w[1]


def test_equal_variance_equal_weights():
    rng = random.Random(3)
    a = [rng.gauss(0, 1) for _ in range(2000)]
    b = [rng.gauss(0, 1) for _ in range(2000)]
    w = inverse_mse_weights([a, b])
    assert abs(w[0] - 0.5) < 0.05


def test_optimal_combination_beats_individual():
    rng = random.Random(1)
    n = 1000
    e1 = [rng.gauss(0, 1) for _ in range(n)]
    e2 = [-0.7 * e1[i] + rng.gauss(0, 0.7) for i in range(n)]   # hedging errors
    w = optimal_combination_weights([e1, e2])
    combo = [w[0] * e1[i] + w[1] * e2[i] for i in range(n)]
    assert _mse(combo) <= min(_mse(e1), _mse(e2)) + 1e-9


def test_combine_applies_weights():
    out = combine_forecasts([[1.0, 1.0], [3.0, 3.0]], [0.25, 0.75])
    assert out == [2.5, 2.5]


def test_validation():
    with pytest.raises(ValueError):
        simple_average_forecast([])
    with pytest.raises(ValueError):
        inverse_mse_weights([])
    with pytest.raises(ValueError):
        optimal_combination_weights([[1.0]])          # < 2 points
    with pytest.raises(ValueError):
        combine_forecasts([[1.0, 2.0]], [0.5, 0.5])   # weight/series mismatch
