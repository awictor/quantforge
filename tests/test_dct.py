"""Discrete cosine transform (DCT-II / DCT-III)."""

import math
import random

import pytest

from quantforge import dct, idct


def test_roundtrip():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(16)]
    xr = idct(dct(x))
    assert max(abs(xr[i] - x[i]) for i in range(16)) < 1e-10


def test_parseval():
    rng = random.Random(5)
    x = [rng.gauss(0, 1) for _ in range(20)]
    X = dct(x)
    assert abs(sum(v * v for v in x) - sum(v * v for v in X)) < 1e-9


def test_constant_only_dc():
    X = dct([5.0] * 8)
    assert abs(X[0] - 40 / math.sqrt(8)) < 1e-9
    assert all(abs(X[k]) < 1e-10 for k in range(1, 8))


def test_energy_compaction_on_smooth():
    smooth = [math.sin(math.pi * i / 32) for i in range(32)]
    X = dct(smooth)
    total = sum(v * v for v in X)
    assert sum(X[k] ** 2 for k in range(4)) / total > 0.95


def test_cosine_mode_isolated():
    n = 16
    mode = [math.cos(math.pi * (2 * i + 1) * 3 / (2 * n)) for i in range(n)]
    X = dct(mode)
    assert max(range(n), key=lambda k: abs(X[k])) == 3


def test_linearity():
    rng = random.Random(7)
    a = [rng.gauss(0, 1) for _ in range(8)]
    b = [rng.gauss(0, 1) for _ in range(8)]
    Xa, Xb = dct(a), dct(b)
    Xab = dct([a[i] + b[i] for i in range(8)])
    assert max(abs(Xab[i] - (Xa[i] + Xb[i])) for i in range(8)) < 1e-12


def test_validation():
    with pytest.raises(ValueError):
        dct([])
    with pytest.raises(ValueError):
        idct([])
