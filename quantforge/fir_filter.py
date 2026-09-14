"""Windowed-sinc FIR filter design and application.

A finite impulse response (FIR) filter with a linear phase is designed by truncating
the ideal filter's (infinite) sinc impulse response to ``numtaps`` and tapering it with
a window to control ripple. This builds lowpass, highpass, and bandpass tap sets and
applies them by convolution. Cutoff frequencies are normalized to the Nyquist rate
(1.0 = half the sampling rate). Pure standard library; uses the library's windows.
"""

import math

from .windows import hann, hamming, blackman


def _window(name, n):
    table = {"hann": hann, "hamming": hamming, "blackman": blackman}
    if name not in table:
        raise ValueError("window must be hann/hamming/blackman")
    return table[name](n)


def _sinc(x):
    if x == 0.0:
        return 1.0
    px = math.pi * x
    return math.sin(px) / px


def fir_lowpass(numtaps, cutoff, window="hamming"):
    """Windowed-sinc lowpass FIR taps.

    ``cutoff`` is the normalized cutoff (0 to 1, fraction of Nyquist). ``numtaps`` should
    be odd for a symmetric (linear-phase, zero-delay-at-center) filter. Returns the tap
    list, normalized to unit DC gain.
    """
    if numtaps < 1:
        raise ValueError("numtaps must be >= 1")
    if not (0.0 < cutoff < 1.0):
        raise ValueError("cutoff must be in (0, 1)")
    w = _window(window, numtaps)
    m = (numtaps - 1) / 2.0
    taps = [cutoff * _sinc(cutoff * (k - m)) * w[k] for k in range(numtaps)]
    s = sum(taps)
    return [t / s for t in taps]


def fir_highpass(numtaps, cutoff, window="hamming"):
    """Windowed-sinc highpass FIR taps (spectral inversion of a lowpass).

    ``numtaps`` must be odd. Unit gain at Nyquist, zero at DC.
    """
    if numtaps % 2 == 0:
        raise ValueError("highpass needs an odd numtaps")
    lp = fir_lowpass(numtaps, cutoff, window)
    mid = (numtaps - 1) // 2
    taps = [-t for t in lp]
    taps[mid] += 1.0
    return taps


def fir_bandpass(numtaps, low, high, window="hamming"):
    """Windowed-sinc bandpass FIR taps passing ``[low, high]`` (normalized cutoffs).

    ``numtaps`` must be odd; the band-pass is the difference of two low-pass filters.
    """
    if numtaps % 2 == 0:
        raise ValueError("bandpass needs an odd numtaps")
    if not (0.0 < low < high < 1.0):
        raise ValueError("need 0 < low < high < 1")
    lp_hi = fir_lowpass(numtaps, high, window)
    lp_lo = fir_lowpass(numtaps, low, window)
    return [lp_hi[i] - lp_lo[i] for i in range(numtaps)]


def fir_apply(taps, x):
    """Apply an FIR filter (tap list) to signal ``x`` by direct convolution.

    Returns the filtered signal, same length as ``x`` (each output is the tap-weighted
    sum of the current and preceding samples; the first ``len(taps)-1`` outputs are the
    transient start-up).
    """
    n = len(x)
    m = len(taps)
    out = []
    for i in range(n):
        acc = 0.0
        for j in range(m):
            if i - j >= 0:
                acc += taps[j] * x[i - j]
        out.append(acc)
    return out
