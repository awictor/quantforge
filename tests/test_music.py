import math
import random

import pytest

from quantforge import music_frequencies, music_pseudospectrum


def close(a, b, tol=0.02):
    return abs(a - b) <= tol


def test_single_tone():
    y = [math.cos(0.8 * k) for k in range(64)]
    f = music_frequencies(y, 2)
    assert len(f) == 1
    assert close(f[0], 0.8, 0.01)


def test_two_separated_tones():
    y = [math.cos(0.5 * k) + 0.8 * math.cos(1.5 * k) for k in range(64)]
    f = music_frequencies(y, 4)
    assert len(f) == 2
    assert close(f[0], 0.5) and close(f[1], 1.5)


def test_closely_spaced_below_fourier_limit():
    # frequency separation 0.06 < Fourier limit 2*pi/64 ~ 0.098
    w1, w2 = 0.60, 0.66
    y = [math.cos(w1 * k) + math.cos(w2 * k) for k in range(64)]
    f = music_frequencies(y, 4, grid=8000)
    assert len(f) == 2
    assert close(f[0], w1) and close(f[1], w2)


def test_noisy_tone():
    random.seed(0)
    y = [math.cos(0.9 * k) + random.uniform(-0.1, 0.1) for k in range(64)]
    assert close(music_frequencies(y, 2)[0], 0.9)


def test_pseudospectrum_peaks_at_tone():
    freqs = [math.pi * i / 2000 for i in range(2001)]
    ps = music_pseudospectrum([math.cos(0.8 * k) for k in range(64)], 2, freqs)
    peak_i = max(range(len(ps)), key=lambda i: ps[i])
    assert close(freqs[peak_i], 0.8, 0.01)


def test_errors():
    with pytest.raises(ValueError):
        music_frequencies([math.cos(k) for k in range(64)], 80)
