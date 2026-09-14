"""Total-variation denoising (Condat's 1-D algorithm)."""

import random

import pytest

from quantforge import tv_denoise, tv_total_variation


def _obj(x, y, lam):
    return (0.5 * sum((x[i] - y[i]) ** 2 for i in range(len(y)))
            + lam * tv_total_variation(x))


def test_lambda_zero_is_identity():
    y = [1.0, 5.0, 2.0, 8.0]
    assert tv_denoise(y, 0.0) == y


def test_large_lambda_collapses_to_mean():
    y = [1.0, 5.0, 2.0, 8.0]
    d = tv_denoise(y, 1e6)
    assert max(d) - min(d) < 1e-6
    assert abs(d[0] - sum(y) / len(y)) < 1e-4


def test_step_is_preserved():
    rng = random.Random(3)
    step = [0.0] * 30 + [10.0] * 30
    noisy = [v + rng.gauss(0.0, 0.5) for v in step]
    d = tv_denoise(noisy, 3.0)
    assert abs(sum(d[:25]) / 25) < 0.5           # low plateau near 0
    assert abs(sum(d[35:]) / 25 - 10.0) < 0.5    # high plateau near 10
    # denoised total variation is far below the noisy input's
    assert tv_total_variation(d) < 0.3 * tv_total_variation(noisy)


def test_result_is_the_optimizer():
    rng = random.Random(7)
    y = [rng.gauss(0.0, 1.0) for _ in range(40)]
    lam = 2.0
    d = tv_denoise(y, lam)
    base = _obj(d, y, lam)
    # no local coordinate perturbation lowers the objective
    for _ in range(3000):
        x = list(d)
        i = rng.randrange(len(x))
        x[i] += rng.uniform(-0.5, 0.5)
        assert _obj(x, y, lam) >= base - 1e-9


def test_matches_independent_iterative_solver():
    y = [1.0, 1.2, 0.9, 5.0, 5.1, 4.8]
    lam = 0.5
    d = tv_denoise(y, lam)
    # a slow subgradient descent should not beat the exact solver's objective
    x = list(y)
    for _ in range(200000):
        g = [x[i] - y[i] for i in range(len(y))]
        for k in range(len(y) - 1):
            g[k] += lam * (1.0 if x[k] > x[k + 1] else (-1.0 if x[k] < x[k + 1] else 0.0))
            g[k + 1] += lam * (1.0 if x[k + 1] > x[k] else (-1.0 if x[k + 1] < x[k] else 0.0))
        for i in range(len(y)):
            x[i] -= 1e-3 * g[i]
    assert _obj(d, y, lam) <= _obj(x, y, lam) + 1e-3


def test_single_element():
    assert tv_denoise([4.0], 5.0) == [4.0]


def test_validation():
    with pytest.raises(ValueError):
        tv_denoise([], 1.0)
    with pytest.raises(ValueError):
        tv_denoise([1.0, 2.0], -1.0)
