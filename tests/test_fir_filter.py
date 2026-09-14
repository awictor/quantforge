"""Windowed-sinc FIR filter design and application."""

import math

import pytest

from quantforge import fir_lowpass, fir_highpass, fir_bandpass, fir_apply


def _gain(taps, freq, N=800):
    y = fir_apply(taps, [math.cos(2 * math.pi * freq * n) for n in range(N)])
    return max(abs(v) for v in y[-100:])


def test_lowpass_unit_dc_gain_and_symmetric():
    lp = fir_lowpass(51, 0.2)
    assert abs(sum(lp) - 1.0) < 1e-9
    assert all(abs(lp[i] - lp[-1 - i]) < 1e-12 for i in range(51))


def test_lowpass_passes_low_blocks_high():
    lp = fir_lowpass(51, 0.2)     # cutoff 0.2 Nyquist = 0.1 cyc/sample
    assert _gain(lp, 0.025) > 0.95      # well below cutoff
    assert _gain(lp, 0.2) < 0.05        # well above


def test_highpass_blocks_dc():
    hp = fir_highpass(51, 0.2)
    assert abs(sum(hp)) < 1e-9          # zero DC gain
    assert _gain(hp, 0.025) < 0.05      # blocks low
    assert _gain(hp, 0.225) > 0.9       # passes high


def test_bandpass_passes_middle():
    bp = fir_bandpass(101, 0.15, 0.35)  # [0.075, 0.175] cyc/sample
    assert _gain(bp, 0.03) < 0.05       # below band
    assert _gain(bp, 0.125) > 0.9       # in band
    assert _gain(bp, 0.25) < 0.05       # above band


def test_two_tone_lowpass_removes_high():
    N = 500
    sig = [math.cos(2 * math.pi * 0.025 * n) + math.cos(2 * math.pi * 0.2 * n)
           for n in range(N)]
    filt = fir_apply(fir_lowpass(101, 0.1), sig)
    assert max(abs(filt[i]) for i in range(400, 500)) < 1.1   # ~ just the low tone


def test_validation():
    with pytest.raises(ValueError):
        fir_lowpass(51, 1.5)
    with pytest.raises(ValueError):
        fir_highpass(50, 0.2)               # even numtaps
    with pytest.raises(ValueError):
        fir_bandpass(101, 0.4, 0.2)         # low > high
