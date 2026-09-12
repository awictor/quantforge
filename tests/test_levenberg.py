"""Levenberg-Marquardt nonlinear least squares."""

import math
import random

import pytest

from quantforge import levenberg_marquardt as lm


def _exp_model(beta, x):
    return beta[0] * math.exp(beta[1] * x)


def test_recovers_exponential_parameters():
    xs = [i * 0.2 for i in range(30)]
    ys = [2.0 * math.exp(0.5 * x) for x in xs]
    r = lm(_exp_model, xs, ys, [1.0, 0.1])
    a, b = r["parameters"]
    assert abs(a - 2.0) < 1e-3 and abs(b - 0.5) < 1e-3
    assert r["residual"] < 1e-6
    assert r["converged"]


def test_noisy_fit_close():
    rng = random.Random(2)
    xs = [i * 0.2 for i in range(30)]
    ys = [2.0 * math.exp(0.5 * x) + rng.gauss(0, 0.05) for x in xs]
    a, b = lm(_exp_model, xs, ys, [1.0, 0.1])["parameters"]
    assert abs(a - 2.0) < 0.1 and abs(b - 0.5) < 0.05


def test_converges_from_bad_start():
    xs = [i * 0.2 for i in range(30)]
    ys = [2.0 * math.exp(0.5 * x) for x in xs]
    a, b = lm(_exp_model, xs, ys, [10.0, -1.0])["parameters"]
    assert abs(a - 2.0) < 0.01 and abs(b - 0.5) < 0.01


def test_linear_model_exact():
    lin = lambda beta, x: beta[0] * x + beta[1]
    r = lm(lin, [0, 1, 2, 3, 4], [1, 4, 7, 10, 13], [0.0, 0.0])
    assert abs(r["parameters"][0] - 3) < 1e-6
    assert abs(r["parameters"][1] - 1) < 1e-6


def test_validation():
    with pytest.raises(ValueError):
        lm(_exp_model, [1, 2], [1], [1, 1])       # length mismatch
    with pytest.raises(ValueError):
        lm(_exp_model, [1], [1], [1, 1, 1])       # n < p
