import math

import pytest

from quantforge import expected_improvement, bayesian_optimize


def test_expected_improvement_properties():
    assert expected_improvement(1.0, 0.0, 2.0) == 0.0        # no variance -> no EI
    assert expected_improvement(1.0, 1.0, 2.0) > 0           # mean below best, uncertain
    assert expected_improvement(2.0, 4.0, 2.0) > expected_improvement(2.0, 1.0, 2.0)  # more var
    big = expected_improvement(0.0, 0.01, 5.0)
    assert abs(big - 5.0) < 0.2                              # mean far below best -> ~gap


def test_optimize_quadratic():
    cand = [i * 0.1 for i in range(-50, 51)]
    res = bayesian_optimize(lambda x: (x - 2.0) ** 2, cand, n_init=3, n_iter=20,
                            length_scale=1.0, noise=1e-6)
    assert abs(res["best_x"] - 2.0) < 0.15
    assert res["best_y"] < 0.05
    assert res["n_eval"] <= 25          # far fewer than 101 candidates


def test_optimize_multimodal():
    cand = [i * 0.05 for i in range(0, 201)]
    g = lambda x: math.sin(x) + 0.1 * x
    res = bayesian_optimize(g, cand, n_init=5, n_iter=30, length_scale=1.0, noise=1e-6)
    true_min = min(g(x) for x in cand)
    assert res["best_y"] <= true_min + 0.1


def test_optimize_2d():
    cand = [[i * 0.2 - 3, j * 0.2 - 3] for i in range(31) for j in range(31)]
    res = bayesian_optimize(lambda p: (p[0] - 1) ** 2 + (p[1] + 1) ** 2, cand,
                            n_init=4, n_iter=40, length_scale=1.0, noise=1e-6)
    assert res["best_y"] < 0.3


def test_empty_candidates_raises():
    with pytest.raises(ValueError):
        bayesian_optimize(lambda x: x, [])
