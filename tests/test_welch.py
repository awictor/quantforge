"""Welch's power-spectral-density estimate."""

import math
import random
import statistics

import pytest

from quantforge import welch_psd, periodogram


def test_peaks_at_sinusoid_frequency():
    n, f0 = 1024, 0.1
    sig = [math.sin(2 * math.pi * f0 * t) for t in range(n)]
    freqs, power = welch_psd(sig, segment_length=128)
    peak = max(range(len(power)), key=lambda k: power[k])
    assert abs(freqs[peak] - 0.1) < 0.02


def test_lower_variance_than_raw_periodogram():
    rng = random.Random(1)
    wn = [rng.gauss(0, 1) for _ in range(1024)]
    _, wp = welch_psd(wn, segment_length=128)
    _, rp = periodogram(wn)
    assert statistics.pvariance(wp) < statistics.pvariance(rp[1:])


def test_noisy_sinusoid_peak():
    rng = random.Random(2)
    sig = [math.sin(2 * math.pi * 0.2 * t) + 0.5 * rng.gauss(0, 1)
           for t in range(1024)]
    freqs, power = welch_psd(sig, segment_length=256)
    peak = max(range(len(power)), key=lambda k: power[k])
    assert abs(freqs[peak] - 0.2) < 0.02


def test_frequencies_in_range():
    rng = random.Random(3)
    x = [rng.gauss(0, 1) for _ in range(512)]
    freqs, power = welch_psd(x)
    assert freqs[0] == 0.0 and abs(freqs[-1] - 0.5) < 1e-9
    assert all(p >= 0.0 for p in power)


def test_validation():
    with pytest.raises(ValueError):
        welch_psd([1.0] * 8)               # too short
    with pytest.raises(ValueError):
        welch_psd(list(range(100)), segment_length=200)
    with pytest.raises(ValueError):
        welch_psd(list(range(100)), overlap=1.0)
