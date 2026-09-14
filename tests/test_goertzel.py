"""Goertzel single-frequency DFT."""

import cmath
import math
import random

import pytest

from quantforge import goertzel, goertzel_power, fft


def test_matches_fft_all_bins():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(64)]
    X = fft(x)
    assert max(abs(goertzel(x, k) - X[k]) for k in range(64)) < 1e-9


def test_matches_direct_dft():
    rng = random.Random(5)
    x = [rng.gauss(0, 1) for _ in range(32)]
    dft = lambda k: sum(x[n] * cmath.exp(-2j * cmath.pi * k * n / len(x))
                        for n in range(len(x)))
    for k in (0, 1, 7, 31):
        assert abs(goertzel(x, k) - dft(k)) < 1e-9


def test_power_equals_coefficient_magnitude_squared():
    rng = random.Random(7)
    x = [rng.gauss(0, 1) for _ in range(50)]
    assert abs(goertzel_power(x, 5) - abs(goertzel(x, 5)) ** 2) < 1e-8


def test_tone_peaks_at_its_bin():
    N, k0 = 64, 8
    tone = [math.cos(2 * math.pi * k0 * n / N) for n in range(N)]
    powers = [goertzel_power(tone, k) for k in range(N)]
    peak = max(range(N), key=lambda k: powers[k])
    assert peak in (k0, N - k0)              # real signal -> mirror bins


def test_dc_bin_is_sum():
    rng = random.Random(9)
    x = [rng.gauss(0, 1) for _ in range(40)]
    assert abs(goertzel(x, 0).real - sum(x)) < 1e-9


def test_absent_frequency_low_power():
    N = 64
    tone = [math.cos(2 * math.pi * 8 * n / N) for n in range(N)]
    assert goertzel_power(tone, 20) < 1e-6 * goertzel_power(tone, 8)


def test_validation():
    with pytest.raises(ValueError):
        goertzel([], 0)
    with pytest.raises(ValueError):
        goertzel([1.0, 2.0], 5)
