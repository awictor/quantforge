"""FFT-based convolution and autocorrelation."""

import random

import pytest

from quantforge import convolve, fft_autocorrelation


def _direct_conv(a, b):
    out = [0.0] * (len(a) + len(b) - 1)
    for i in range(len(a)):
        for j in range(len(b)):
            out[i + j] += a[i] * b[j]
    return out


def test_matches_direct_convolution():
    a, b = [1, 2, 3], [4, 5, 6, 7]
    c = convolve(a, b)
    d = _direct_conv(a, b)
    assert max(abs(c[i] - d[i]) for i in range(len(c))) < 1e-9


def test_polynomial_product():
    # (1 + x)(1 + x) = 1 + 2x + x^2
    c = convolve([1, 1], [1, 1])
    assert [round(v, 6) for v in c] == [1.0, 2.0, 1.0]


def test_delta_identity():
    c = convolve([1, 0, 0], [5, 6, 7])
    assert [round(c[i], 6) for i in range(3)] == [5.0, 6.0, 7.0]


def test_output_length():
    assert len(convolve([1, 2], [3, 4, 5])) == 4


def test_autocorrelation_matches_direct():
    rng = random.Random(1)
    x = [rng.gauss(0, 1) for _ in range(200)]
    acf = fft_autocorrelation(x, max_lag=5)

    def direct(lag):
        n = len(x)
        m = sum(x) / n
        c0 = sum((v - m) ** 2 for v in x)
        ck = sum((x[i] - m) * (x[i - lag] - m) for i in range(lag, n))
        return ck / c0

    assert abs(acf[0] - 1.0) < 1e-9
    assert max(abs(acf[k] - direct(k)) for k in range(6)) < 1e-9


def test_validation():
    with pytest.raises(ValueError):
        convolve([], [1.0])
    with pytest.raises(ValueError):
        fft_autocorrelation([1.0])
    with pytest.raises(ValueError):
        fft_autocorrelation([2.0, 2.0, 2.0])   # zero variance
