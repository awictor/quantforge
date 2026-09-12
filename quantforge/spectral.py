"""Spectral analysis of a real signal: DFT, periodogram, dominant frequency.

Cyclical structure in a series -- a seasonal pattern, a dominant trading cycle --
shows up as concentrated power at particular frequencies. The discrete Fourier
transform decomposes the signal into sinusoids; the periodogram is the squared
magnitude of that transform and estimates how variance is distributed across
frequency. A direct O(n^2) DFT is used (no radix-2 restriction), which is plenty
for the series lengths this is used on. Pure standard library.
"""

import cmath
import math


def dft(x):
    """Discrete Fourier transform of a real (or complex) sequence.

    Returns the ``n`` complex coefficients ``X_k = sum_t x_t exp(-2 pi i k t / n)``.
    """
    n = len(x)
    if n == 0:
        raise ValueError("empty sequence")
    out = []
    for k in range(n):
        acc = 0j
        ang = -2.0 * math.pi * k / n
        for t in range(n):
            acc += x[t] * cmath.exp(1j * ang * t)
        out.append(acc)
    return out


def periodogram(x):
    """One-sided periodogram of a real series.

    Returns ``(freqs, power)`` where ``freqs`` are normalized frequencies in
    cycles per sample over ``[0, 0.5]`` and ``power[k] = |X_k|^2 / n``. The DC
    term is included at frequency 0. Peaks mark dominant cycles.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least 2 points")
    coeffs = dft(x)
    half = n // 2
    freqs = []
    power = []
    for k in range(half + 1):
        freqs.append(k / n)
        power.append(abs(coeffs[k]) ** 2 / n)
    return freqs, power


def dominant_frequency(x):
    """Frequency (cycles per sample) of the largest non-DC periodogram peak.

    Ignores the zero-frequency (mean) component. Returns the frequency; its
    reciprocal is the dominant period in samples.
    """
    freqs, power = periodogram(x)
    best_k = None
    best_p = -1.0
    for k in range(1, len(power)):        # skip DC
        if power[k] > best_p:
            best_p = power[k]
            best_k = k
    return freqs[best_k]


def spectral_energy(x):
    """Total periodogram energy, ``sum_k |X_k|^2 / n`` over all n frequencies.

    By Parseval's theorem this equals ``sum_t x_t^2`` (the time-domain energy),
    which the tests use as a consistency check.
    """
    n = len(x)
    coeffs = dft(x)
    return sum(abs(c) ** 2 for c in coeffs) / n
