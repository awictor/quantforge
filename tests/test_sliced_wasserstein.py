import math
import random

import pytest

from quantforge import sliced_wasserstein
from quantforge.wasserstein import wasserstein_distance


def close(a, b, tol=1e-2):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def test_one_dimensional_reduces_to_wasserstein():
    random.seed(0)
    x = [random.gauss(0, 1) for _ in range(200)]
    y = [random.gauss(1, 1) for _ in range(200)]
    sw = sliced_wasserstein(x, y, n_projections=300, p=2, seed=1)
    assert close(sw, wasserstein_distance(x, y, p=2), 1e-6)


def _shift_setup(d=3, n=400, t=(2.0, -1.0, 0.5), seed=2):
    random.seed(seed)
    base = [[random.gauss(0, 1) for _ in range(d)] for _ in range(n)]
    shifted = [[base[i][k] + t[k] for k in range(d)] for i in range(n)]
    return base, shifted, t, d


def test_translation_matches_analytic():
    base, shifted, t, d = _shift_setup()
    sw = sliced_wasserstein(base, shifted, n_projections=2000, p=2, seed=3)
    expect = math.sqrt(sum(ti * ti for ti in t)) / math.sqrt(d)
    assert close(sw, expect, 5e-2)


def test_identical_is_zero():
    base, _, _, _ = _shift_setup()
    assert sliced_wasserstein(base, base, n_projections=100, seed=5) < 1e-9


def test_reproducible_with_seed():
    base, shifted, _, _ = _shift_setup()
    a = sliced_wasserstein(base, shifted, n_projections=200, seed=7)
    b = sliced_wasserstein(base, shifted, n_projections=200, seed=7)
    assert a == b


def test_symmetry():
    base, shifted, _, _ = _shift_setup()
    a = sliced_wasserstein(base, shifted, n_projections=200, seed=9)
    b = sliced_wasserstein(shifted, base, n_projections=200, seed=9)
    assert close(a, b, 1e-12)


def test_triangle_inequality():
    base, shifted, t, d = _shift_setup()
    n = len(base)
    z = [[base[i][k] + 5 * t[k] for k in range(d)] for i in range(n)]
    sxz = sliced_wasserstein(base, z, 3000, seed=11)
    sxy = sliced_wasserstein(base, shifted, 3000, seed=11)
    syz = sliced_wasserstein(shifted, z, 3000, seed=11)
    assert sxz <= sxy + syz + 1e-6


def test_empty_raises():
    with pytest.raises(ValueError):
        sliced_wasserstein([], [1.0], 10)
