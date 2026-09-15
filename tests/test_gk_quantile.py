import math
import random

import pytest

from quantforge import GKQuantile


def _exact_rank(data, x):
    return sum(1 for v in data if v <= x)


def _data(n=100000, seed=0):
    random.seed(seed)
    return [random.gauss(0, 1) for _ in range(n)]


def test_small_stream_bound():
    eps = 0.02
    data = _data(10000, seed=5)
    n = len(data)
    gk = GKQuantile(eps).add_all(data)
    assert gk.n == n
    for q in (0.25, 0.5, 0.75):
        assert abs(_exact_rank(data, gk.quantile(q)) - q * n) / n <= eps + 2e-3


@pytest.mark.slow
def test_deterministic_rank_guarantee():
    eps = 0.01
    data = _data()
    n = len(data)
    gk = GKQuantile(eps).add_all(data)
    assert gk.n == n
    maxerr = max(abs(_exact_rank(data, gk.quantile(q)) - q * n) / n
                 for q in (0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999))
    assert maxerr <= eps + 1e-3


@pytest.mark.slow
def test_sorted_input_still_bounded():
    eps = 0.01
    data = _data()
    n = len(data)
    gk = GKQuantile(eps).add_all(sorted(data))
    for q in (0.1, 0.5, 0.9):
        assert abs(_exact_rank(data, gk.quantile(q)) - q * n) / n <= eps + 1e-3


@pytest.mark.slow
def test_space_is_bounded():
    eps = 0.01
    data = _data()
    gk = GKQuantile(eps).add_all(data)
    assert len(gk._tuples) < (1 / eps) * math.log(eps * len(data)) * 3


@pytest.mark.slow
def test_uniform_quantiles():
    random.seed(1)
    gk = GKQuantile(0.01).add_all([random.random() for _ in range(100000)])
    for q in (0.1, 0.3, 0.5, 0.7, 0.9):
        assert abs(gk.quantile(q) - q) < 0.02


def test_epsilon_tradeoff():
    data = _data(50000)
    coarse = GKQuantile(0.05).add_all(data)
    fine = GKQuantile(0.005).add_all(data)
    assert len(fine._tuples) > len(coarse._tuples)


def test_errors():
    with pytest.raises(ValueError):
        GKQuantile(0)
    with pytest.raises(ValueError):
        GKQuantile(0.01).quantile(0.5)
