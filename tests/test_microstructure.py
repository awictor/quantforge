"""Order-flow microstructure measures."""

import random

import pytest

from quantforge import kyle_lambda_regression, order_flow_imbalance, vpin


def test_kyle_lambda_recovers_known_impact():
    lam = 0.002
    rng = random.Random(5)
    q = [rng.gauss(0, 1000) for _ in range(3000)]
    dp = [lam * qi + rng.gauss(0, 0.05) for qi in q]
    assert abs(kyle_lambda_regression(dp, q) - lam) < 1e-4


def test_kyle_lambda_zero_when_no_impact():
    rng = random.Random(1)
    q = [rng.gauss(0, 1000) for _ in range(2000)]
    dp = [rng.gauss(0, 0.05) for _ in range(2000)]   # independent of flow
    assert abs(kyle_lambda_regression(dp, q)) < 1e-4


def test_order_flow_imbalance_sign():
    assert order_flow_imbalance([100, 120], [50, 40]) > 0
    assert order_flow_imbalance([40, 50], [100, 120]) < 0
    assert order_flow_imbalance([100, 100], [100, 100]) == 0.0


def test_order_flow_imbalance_bounds():
    assert order_flow_imbalance([100], [0]) == 1.0
    assert order_flow_imbalance([0], [100]) == -1.0


def test_vpin_balanced_is_zero():
    assert vpin([100, 100, 100], [100, 100, 100]) == 0.0


def test_vpin_one_sided_is_one():
    assert vpin([100, 100], [0, 0]) == 1.0


def test_vpin_in_unit_interval():
    v = vpin([80, 60, 90], [20, 40, 10])
    assert 0.0 <= v <= 1.0


def test_vpin_skips_empty_buckets():
    # A zero-volume bucket is ignored, not counted as balanced.
    v = vpin([100, 0], [0, 0])
    assert v == 1.0


def test_validation():
    with pytest.raises(ValueError):
        kyle_lambda_regression([0.1, 0.2], [1.0, 2.0])       # < 3 obs
    with pytest.raises(ValueError):
        kyle_lambda_regression([0.1, 0.2, 0.3], [5.0, 5.0, 5.0])  # constant flow
    with pytest.raises(ValueError):
        order_flow_imbalance([100], [-1])
    with pytest.raises(ValueError):
        vpin([0, 0], [0, 0])                                 # all empty
