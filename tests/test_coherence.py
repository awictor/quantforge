"""Cross-spectral density and magnitude-squared coherence (Welch)."""

import math
import random

import pytest

from quantforge import cross_spectral_density, coherence


def test_identical_signals_full_coherence():
    rng = random.Random(1)
    x = [math.sin(2 * math.pi * 0.1 * n) + rng.gauss(0, 0.5) for n in range(2048)]
    _, c = coherence(x, x, segment_length=256)
    assert min(c) > 0.999


def test_independent_signals_low_coherence():
    rng = random.Random(2)
    a = [rng.gauss(0, 1) for _ in range(2048)]
    b = [rng.gauss(0, 1) for _ in range(2048)]
    _, c = coherence(a, b, segment_length=256)
    assert sum(c) / len(c) < 0.3


def test_linearly_related_high_at_signal_band():
    rng = random.Random(3)
    N = 2048
    base = [math.sin(2 * math.pi * 0.1 * n) for n in range(N)]
    x = [base[n] + rng.gauss(0, 0.2) for n in range(N)]
    y = [2 * base[n] + rng.gauss(0, 0.2) for n in range(N)]
    f, c = coherence(x, y, segment_length=256)
    kbin = min(range(len(f)), key=lambda k: abs(f[k] - 0.1))
    assert c[kbin] > 0.9
    assert all(0.0 <= v <= 1.0 + 1e-9 for v in c)


def test_cross_spectrum_peaks_at_shared_frequency():
    rng = random.Random(4)
    N = 2048
    base = [math.sin(2 * math.pi * 0.1 * n) for n in range(N)]
    x = [base[n] + rng.gauss(0, 0.2) for n in range(N)]
    y = [2 * base[n] + rng.gauss(0, 0.2) for n in range(N)]
    f, xs = cross_spectral_density(x, y, segment_length=256)
    mags = [abs(v) for v in xs]
    peak = max(range(len(mags)), key=lambda k: mags[k])
    assert abs(f[peak] - 0.1) < 0.02


def test_validation():
    with pytest.raises(ValueError):
        coherence([1, 2, 3], [1, 2, 3])
    with pytest.raises(ValueError):
        coherence([1] * 20, [1] * 19)
