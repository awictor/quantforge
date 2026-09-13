"""Diebold-Mariano equal-predictive-accuracy test."""

import random

import pytest

from quantforge import diebold_mariano


def test_favors_more_accurate_forecast():
    rng = random.Random(1)
    e1 = [rng.gauss(0, 0.5) for _ in range(500)]     # accurate
    e2 = [rng.gauss(0, 1.5) for _ in range(500)]     # worse
    dm, p = diebold_mariano(e1, e2)
    assert dm < 0            # forecast 1 has smaller loss
    assert p < 0.001


def test_antisymmetric_in_arguments():
    rng = random.Random(1)
    e1 = [rng.gauss(0, 0.5) for _ in range(500)]
    e2 = [rng.gauss(0, 1.5) for _ in range(500)]
    dm_ab, p_ab = diebold_mariano(e1, e2)
    dm_ba, p_ba = diebold_mariano(e2, e1)
    assert abs(dm_ab + dm_ba) < 1e-9
    assert abs(p_ab - p_ba) < 1e-12


def test_equal_accuracy_not_rejected():
    rng = random.Random(2)
    a = [rng.gauss(0, 1) for _ in range(500)]
    b = [rng.gauss(0, 1) for _ in range(500)]
    _, p = diebold_mariano(a, b)
    assert p > 0.05


def test_absolute_loss_option():
    rng = random.Random(3)
    e1 = [rng.gauss(0, 0.5) for _ in range(400)]
    e2 = [rng.gauss(0, 1.2) for _ in range(400)]
    dm, p = diebold_mariano(e1, e2, power=1)
    assert dm < 0 and p < 0.01


def test_horizon_lags():
    rng = random.Random(4)
    e1 = [rng.gauss(0, 0.5) for _ in range(400)]
    e2 = [rng.gauss(0, 1.2) for _ in range(400)]
    # A larger horizon uses more HAC lags but should still favor forecast 1.
    dm, p = diebold_mariano(e1, e2, h=5)
    assert dm < 0


def test_validation():
    with pytest.raises(ValueError):
        diebold_mariano([0.1], [0.2])
    with pytest.raises(ValueError):
        diebold_mariano([0.1, 0.2], [0.3], h=1)
    with pytest.raises(ValueError):
        diebold_mariano([0.1, 0.2], [0.3, 0.4], h=0)
    with pytest.raises(ValueError):
        diebold_mariano([0.1, 0.2], [0.1, 0.2])       # identical -> zero variance
