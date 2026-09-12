"""Spectral analysis: DFT, periodogram, dominant frequency."""

import math
import random

import pytest

from quantforge import dft, periodogram, dominant_frequency, spectral_energy


def _sine(n, cycles, amp=1.0):
    return [amp * math.sin(2 * math.pi * cycles * t / n) for t in range(n)]


def test_pure_sine_dominant_frequency():
    n, f0 = 128, 8
    x = _sine(n, f0)
    assert abs(dominant_frequency(x) - f0 / n) < 1e-9
    assert abs(1.0 / dominant_frequency(x) - n / f0) < 1e-9


def test_parseval_energy_conservation():
    x = _sine(128, 8)
    time_energy = sum(v * v for v in x)
    assert abs(spectral_energy(x) - time_energy) < 1e-6


def test_constant_signal_is_dc_only():
    freqs, power = periodogram([3.0] * 128)
    assert power[0] > 1e-6
    assert all(p < 1e-9 for p in power[1:])


def test_white_noise_spread_across_frequencies():
    rng = random.Random(1)
    wn = [rng.gauss(0, 1) for _ in range(256)]
    _, power = periodogram(wn)
    total = sum(power[1:])
    assert max(power[1:]) / total < 0.15      # no single frequency dominates


def test_two_tone_picks_stronger():
    n = 128
    y = [3 * math.sin(2 * math.pi * 4 * t / n) + 1 * math.sin(2 * math.pi * 20 * t / n)
         for t in range(n)]
    assert abs(dominant_frequency(y) - 4 / n) < 1e-9


def test_dft_inverts_via_energy():
    x = _sine(64, 5, amp=2.0)
    coeffs = dft(x)
    assert len(coeffs) == len(x)


def test_validation():
    with pytest.raises(ValueError):
        dft([])
    with pytest.raises(ValueError):
        periodogram([1.0])
