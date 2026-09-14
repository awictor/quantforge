"""Butterworth IIR filter design (bilinear transform, biquad cascade)."""

import cmath
import math

import pytest

from quantforge import (
    butter_lowpass,
    butter_highpass,
    sosfilt,
    iir_frequency_response,
)


def _stable(sections):
    for b, a in sections:
        a = (a + [0.0, 0.0, 0.0])[:3]
        if a[2] != 0.0:
            disc = a[1] * a[1] - 4 * a[0] * a[2]
            r1 = (-a[1] + cmath.sqrt(disc)) / (2 * a[0])
            r2 = (-a[1] - cmath.sqrt(disc)) / (2 * a[0])
            if abs(r1) >= 1.0 - 1e-9 or abs(r2) >= 1.0 - 1e-9:
                return False
        elif a[1] != 0.0:
            if abs(-a[2] / a[1]) >= 1.0:
                return False
    return True


def test_lowpass_gain_landmarks():
    lp = butter_lowpass(4, 0.1)
    H = iir_frequency_response(lp, [0.0, 0.1, 0.49])
    assert abs(H[0] - 1.0) < 1e-6            # unit DC gain
    assert abs(H[1] - 0.7071) < 1e-3         # -3 dB at cutoff
    assert H[2] < 0.01                       # deep stopband at Nyquist


def test_lowpass_monotone_rolloff():
    lp = butter_lowpass(4, 0.1)
    grid = [i / 200 for i in range(1, 100)]
    H = iir_frequency_response(lp, grid)
    assert all(H[i] >= H[i + 1] - 1e-9 for i in range(len(H) - 1))


def test_highpass_gain_landmarks():
    hp = butter_highpass(4, 0.1)
    H = iir_frequency_response(hp, [0.0, 0.1, 0.49])
    assert H[0] < 0.01                       # near-zero DC gain
    assert abs(H[1] - 0.7071) < 1e-3         # -3 dB at cutoff
    assert abs(H[2] - 1.0) < 1e-6            # unit gain at Nyquist


def test_time_domain_separation():
    N = 2000
    low = [math.cos(2 * math.pi * 0.02 * n) for n in range(N)]
    high = [math.cos(2 * math.pi * 0.3 * n) for n in range(N)]
    lp = butter_lowpass(4, 0.1)
    assert max(abs(v) for v in sosfilt(lp, low)[500:]) > 0.95    # passes low
    assert max(abs(v) for v in sosfilt(lp, high)[500:]) < 0.05   # blocks high
    sig = [low[n] + high[n] for n in range(N)]
    assert max(abs(v) for v in sosfilt(lp, sig)[500:]) < 1.1     # ~ low tone alone


def test_stability_even_and_odd_order():
    assert _stable(butter_lowpass(4, 0.1))
    assert _stable(butter_highpass(4, 0.1))
    lp3 = butter_lowpass(3, 0.15)
    assert _stable(lp3)
    assert abs(iir_frequency_response(lp3, [0.0])[0] - 1.0) < 1e-6


def test_higher_order_is_sharper():
    grid = [i / 200 for i in range(1, 100)]
    h2 = iir_frequency_response(butter_lowpass(2, 0.1), grid)
    h8 = iir_frequency_response(butter_lowpass(8, 0.1), grid)
    # at a frequency well into the stopband the higher order attenuates more
    idx = grid.index(0.2)
    assert h8[idx] < h2[idx]


def test_validation():
    with pytest.raises(ValueError):
        butter_lowpass(0, 0.1)               # order < 1
    with pytest.raises(ValueError):
        butter_lowpass(2, 0.6)               # cutoff >= 0.5
    with pytest.raises(ValueError):
        butter_highpass(2, 0.0)              # cutoff <= 0
