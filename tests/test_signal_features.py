"""Time- and frequency-domain signal features."""

import math
import random

import pytest

from quantforge import (
    zero_crossing_rate,
    rms,
    crest_factor,
    spectral_centroid,
    spectral_bandwidth,
    spectral_flatness,
)

N = 2048


def _sine(f=0.05):
    return [math.sin(2 * math.pi * f * n) for n in range(N)]


def test_zero_crossing_rate():
    assert abs(zero_crossing_rate([(-1) ** i for i in range(100)]) - 1.0) < 1e-9
    assert zero_crossing_rate([5.0] * 100) == 0.0


def test_rms_and_crest_of_sine():
    s = _sine()
    assert abs(rms(s) - 1 / math.sqrt(2)) < 1e-3
    assert abs(crest_factor(s) - math.sqrt(2)) < 1e-2
    sq = [1.0 if math.sin(2 * math.pi * 0.05 * n) >= 0 else -1.0 for n in range(N)]
    assert abs(crest_factor(sq) - 1.0) < 1e-9


def test_spectral_centroid_of_tone():
    tone = [math.cos(2 * math.pi * 0.2 * n) for n in range(N)]
    assert abs(spectral_centroid(tone) - 0.2) < 0.01


def test_bandwidth_tone_less_than_noise():
    rng = random.Random(1)
    tone = [math.cos(2 * math.pi * 0.2 * n) for n in range(N)]
    noise = [rng.gauss(0, 1) for _ in range(N)]
    assert spectral_bandwidth(tone) < spectral_bandwidth(noise)


def test_flatness_tone_vs_noise():
    rng = random.Random(2)
    tone = [math.cos(2 * math.pi * 0.2 * n) for n in range(N)]
    noise = [rng.gauss(0, 1) for _ in range(N)]
    assert spectral_flatness(tone) < 0.1
    assert spectral_flatness(noise) > 0.3


def test_validation():
    with pytest.raises(ValueError):
        rms([])
    with pytest.raises(ValueError):
        zero_crossing_rate([1])
    with pytest.raises(ValueError):
        crest_factor([0, 0, 0])
