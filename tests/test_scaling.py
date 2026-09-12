"""Feature scaling: standardize, min-max, robust."""

import math
import random

import pytest

from quantforge import (
    fit_standardize, fit_min_max, fit_robust,
    scale_transform, scale_inverse_transform,
)


def _data(seed):
    rng = random.Random(seed)
    return [[rng.gauss(5, 2), rng.gauss(-3, 4)] for _ in range(1000)]


def test_standardize_mean_zero_std_one():
    X = _data(1)
    Z = scale_transform(fit_standardize(X), X)
    col = [z[0] for z in Z]
    m = sum(col) / len(col)
    s = math.sqrt(sum((v - m) ** 2 for v in col) / len(col))
    assert abs(m) < 1e-9
    assert abs(s - 1.0) < 1e-9


def test_standardize_round_trip():
    X = _data(1)
    p = fit_standardize(X)
    back = scale_inverse_transform(p, scale_transform(p, X))
    assert max(abs(back[i][0] - X[i][0]) for i in range(len(X))) < 1e-9


def test_min_max_bounds():
    X = _data(1)
    M = scale_transform(fit_min_max(X), X)
    col = [m[0] for m in M]
    assert abs(min(col)) < 1e-12
    assert abs(max(col) - 1.0) < 1e-12


def test_robust_median_resistant_to_outlier():
    X = _data(1)
    clean = fit_robust(X)["center"][0]
    with_outlier = fit_robust(X + [[1000.0, 1000.0]])["center"][0]
    assert abs(with_outlier - clean) < 0.5


def test_transform_uses_fitted_params_on_new_data():
    X = _data(1)
    p = fit_standardize(X)
    # The train mean maps to ~0 under its own params.
    zt = scale_transform(p, [p["center"]])
    assert all(abs(v) < 1e-12 for v in zt[0])


def test_validation():
    p = fit_standardize(_data(1))
    with pytest.raises(ValueError):
        scale_transform(p, [[1.0, 2.0, 3.0]])   # column mismatch
    with pytest.raises(ValueError):
        fit_standardize([])                      # empty
