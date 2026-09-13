"""Radix-2 FFT / inverse FFT."""

import cmath
import math
import random

import pytest

from quantforge import fft, ifft


def _dft(x):
    n = len(x)
    return [sum(x[t] * cmath.exp(-2j * math.pi * k * t / n) for t in range(n))
            for k in range(n)]


def test_matches_direct_dft():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(8)]
    X, Xd = fft(x), _dft(x)
    assert max(abs(X[k] - Xd[k]) for k in range(8)) < 1e-9


def test_roundtrip():
    rng = random.Random(2)
    x = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(16)]
    xr = ifft(fft(x))
    assert max(abs(xr[k] - x[k]) for k in range(16)) < 1e-9


def test_delta_flat_spectrum():
    delta = [1.0] + [0.0] * 7
    assert all(abs(v - 1.0) < 1e-12 for v in fft(delta))


def test_constant_spikes_at_zero():
    Xc = fft([3.0] * 8)
    assert abs(Xc[0].real - 24.0) < 1e-9
    assert all(abs(Xc[k]) < 1e-9 for k in range(1, 8))


def test_linearity():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(8)]
    y = [rng.gauss(0, 1) for _ in range(8)]
    lhs = fft([2 * x[i] + 3 * y[i] for i in range(8)])
    Xf, Yf = fft(x), fft(y)
    rhs = [2 * Xf[i] + 3 * Yf[i] for i in range(8)]
    assert max(abs(lhs[i] - rhs[i]) for i in range(8)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        fft([1.0, 2.0, 3.0])       # length not a power of two
    with pytest.raises(ValueError):
        ifft([])
