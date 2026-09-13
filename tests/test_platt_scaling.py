"""Platt scaling sigmoid calibration."""

import math
import random

import pytest

from quantforge import platt_fit, platt_predict, platt_calibrate


def _smoothed_targets(labels):
    n = len(labels)
    p1 = sum(1 for l in labels if l == 1)
    p0 = n - p1
    hi = (p1 + 1.0) / (p1 + 2.0)
    lo = 1.0 / (p0 + 2.0)
    return [hi if l == 1 else lo for l in labels]


def _sig(z):
    return 1.0 / (1.0 + math.exp(-z)) if z >= 0 else math.exp(z) / (1.0 + math.exp(z))


def _loss(A, B, scores, t):
    L = 0.0
    for i in range(len(scores)):
        p = min(max(_sig(A * scores[i] + B), 1e-15), 1 - 1e-15)
        L -= t[i] * math.log(p) + (1 - t[i]) * math.log(1 - p)
    return L


def test_fit_is_loss_minimum():
    rng = random.Random(2)
    scores = [rng.gauss(0, 1) for _ in range(200)]
    labels = [1 if s + rng.gauss(0, 0.5) > 0 else 0 for s in scores]
    A, B = platt_fit(scores, labels)
    t = _smoothed_targets(labels)
    L0 = _loss(A, B, scores, t)
    # No small perturbation lowers the loss.
    for da in (-0.2, -0.05, 0.05, 0.2):
        for db in (-0.2, -0.05, 0.05, 0.2):
            assert _loss(A + da, B + db, scores, t) >= L0 - 1e-9


def test_gradient_zero_at_fit():
    rng = random.Random(5)
    scores = [rng.gauss(0, 1) for _ in range(150)]
    labels = [1 if s > 0 else 0 for s in scores]
    A, B = platt_fit(scores, labels)
    t = _smoothed_targets(labels)
    gA = sum((_sig(A * scores[i] + B) - t[i]) * scores[i] for i in range(len(scores)))
    gB = sum(_sig(A * scores[i] + B) - t[i] for i in range(len(scores)))
    assert abs(gA) < 1e-6 and abs(gB) < 1e-6


def test_probability_rises_with_score():
    rng = random.Random(7)
    scores = [rng.gauss(0, 1) for _ in range(300)]
    labels = [1 if s > 0 else 0 for s in scores]
    A, B = platt_fit(scores, labels)
    assert A > 0
    probs = platt_predict(sorted(scores), A, B)
    assert all(probs[i] <= probs[i + 1] + 1e-12 for i in range(len(probs) - 1))


def test_predict_in_unit_interval():
    A, B = platt_fit([-2, -1, 0, 1, 2], [0, 0, 0, 1, 1])
    for p in platt_predict([-5, 0, 5], A, B):
        assert 0.0 <= p <= 1.0


def test_calibrate_helper_predict_and_direct():
    rng = random.Random(3)
    scores = [rng.gauss(0, 1) for _ in range(100)]
    labels = [1 if s > 0 else 0 for s in scores]
    (A, B), predict = platt_calibrate(scores, labels)
    (A2, B2), new_probs = platt_calibrate(scores, labels, new_scores=[0.0, 1.0])
    assert (A, B) == (A2, B2)
    assert predict([0.0, 1.0]) == new_probs


def test_validation():
    with pytest.raises(ValueError):
        platt_fit([1, 2], [1])
    with pytest.raises(ValueError):
        platt_fit([1], [2])
