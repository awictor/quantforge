"""Spectral window functions."""

import math

import pytest

from quantforge import (hann, hamming, blackman, bartlett, rectangular, apply_window, fft)


def test_hann_ends_and_center():
    h = hann(9)
    assert h[0] == 0.0 and h[-1] == 0.0
    assert abs(h[4] - 1.0) < 1e-12


def test_hamming_ends():
    hm = hamming(9)
    assert abs(hm[0] - 0.08) < 1e-9 and abs(hm[-1] - 0.08) < 1e-9


def test_bartlett_triangular():
    bt = bartlett(9)
    assert bt[0] == 0.0 and abs(bt[4] - 1.0) < 1e-12 and abs(bt[2] - 0.5) < 1e-12


def test_rectangular_all_ones():
    assert rectangular(4) == [1.0, 1.0, 1.0, 1.0]


def test_all_symmetric():
    for fn in (hann, hamming, blackman, bartlett):
        w = fn(21)
        assert all(abs(w[i] - w[-1 - i]) < 1e-12 for i in range(21))


def test_window_reduces_leakage():
    N = 64
    sig = [math.cos(2 * math.pi * 8.5 * n / N) for n in range(N)]   # off-bin -> leaks
    Xrect = fft(sig)
    Xhann = fft(apply_window(sig, "hann"))
    assert abs(Xhann[30]) < abs(Xrect[30])       # far sidelobe suppressed


def test_apply_custom_window():
    assert apply_window([1, 2, 3], [1, 0, 1]) == [1, 0, 3]


def test_validation():
    with pytest.raises(ValueError):
        hann(0)
    with pytest.raises(ValueError):
        apply_window([1, 2], "unknown")
    with pytest.raises(ValueError):
        apply_window([1, 2, 3], [1, 1])
