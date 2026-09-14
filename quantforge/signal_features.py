"""Time- and frequency-domain signal features.

Compact descriptors that summarize a signal for classification, monitoring, or audio
analysis: the zero-crossing rate and RMS/crest factor in the time domain, and the
spectral centroid, bandwidth, and flatness in the frequency domain. The spectral features
are computed from the periodogram (magnitude spectrum). Pure standard library on the DFT.
"""

import math

from .spectral import dft


def zero_crossing_rate(x):
    """Fraction of adjacent sample pairs that straddle zero (sign changes / (n-1)).

    A rough pitch/noisiness proxy: high for noisy or high-frequency signals, low for
    smooth low-frequency ones. Returns a value in ``[0, 1]``.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least two samples")
    crossings = 0
    for i in range(1, n):
        if (x[i - 1] < 0) != (x[i] < 0):
            crossings += 1
    return crossings / (n - 1)


def rms(x):
    """Root-mean-square amplitude of a signal."""
    if not x:
        raise ValueError("input must be non-empty")
    return math.sqrt(sum(v * v for v in x) / len(x))


def crest_factor(x):
    """Crest factor: peak amplitude divided by RMS.

    High for impulsive/peaky signals (a lone spike), low for signals that fill their
    range (~1.41 for a sine, 1.0 for a square wave). Zero-signal raises.
    """
    r = rms(x)
    if r == 0.0:
        raise ValueError("signal is all zeros")
    return max(abs(v) for v in x) / r


def _power_spectrum(x):
    """One-sided power spectrum and its bin frequencies (cycles/sample)."""
    n = len(x)
    if n < 2:
        raise ValueError("need at least two samples")
    m = sum(x) / n
    coeffs = dft([v - m for v in x])
    half = n // 2
    power = [abs(coeffs[k]) ** 2 for k in range(half + 1)]
    freqs = [k / n for k in range(half + 1)]
    return freqs, power


def spectral_centroid(x):
    """Spectral centroid: the power-weighted mean frequency (cycles/sample).

    The spectrum's "center of mass" -- higher for brighter/higher-pitched signals. Zero
    for a DC-only (constant) signal. Computed from the mean-removed periodogram.
    """
    freqs, power = _power_spectrum(x)
    total = sum(power)
    if total == 0.0:
        return 0.0
    return sum(freqs[k] * power[k] for k in range(len(power))) / total


def spectral_bandwidth(x):
    """Spectral bandwidth: the power-weighted standard deviation about the centroid.

    Measures how spread out the spectrum is (narrow for a pure tone, wide for noise).
    """
    freqs, power = _power_spectrum(x)
    total = sum(power)
    if total == 0.0:
        return 0.0
    c = sum(freqs[k] * power[k] for k in range(len(power))) / total
    var = sum(power[k] * (freqs[k] - c) ** 2 for k in range(len(power))) / total
    return math.sqrt(var)


def spectral_flatness(x):
    """Spectral flatness (Wiener entropy): geometric mean / arithmetic mean of the spectrum.

    Near ``1`` for white-noise-like flat spectra, near ``0`` for tonal signals with power
    concentrated in a few bins. Uses the mean-removed positive-frequency power bins.
    """
    freqs, power = _power_spectrum(x)
    p = [v for v in power[1:] if v > 0]        # skip DC; drop exact zeros for the geo-mean
    if not p:
        return 0.0
    geo = math.exp(sum(math.log(v) for v in p) / len(p))
    arith = sum(p) / len(p)
    return geo / arith if arith > 0 else 0.0
