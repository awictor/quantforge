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


def welch_psd(x, segment_length=None, overlap=0.5):
    """Welch's power-spectral-density estimate: averaged windowed periodograms.

    Splits ``x`` into overlapping segments of ``segment_length`` (default ``n // 8``,
    clamped to at least 8), applies a Hann window to each, and averages their
    periodograms. Averaging trades frequency resolution for a much lower-variance
    spectral estimate than the raw periodogram. ``overlap`` is the fractional segment
    overlap in ``[0, 1)``. Returns ``(freqs, power)`` with one-sided normalized
    frequencies in ``[0, 0.5]``.
    """
    n = len(x)
    if n < 16:
        raise ValueError("need at least 16 points")
    if segment_length is None:
        segment_length = max(8, n // 8)
    if segment_length > n:
        raise ValueError("segment_length exceeds the series length")
    if not (0.0 <= overlap < 1.0):
        raise ValueError("overlap must be in [0, 1)")

    L = segment_length
    step = max(1, int(L * (1.0 - overlap)))
    # Hann window and its power (for normalization).
    win = [0.5 - 0.5 * math.cos(2.0 * math.pi * i / (L - 1)) for i in range(L)]
    win_power = sum(w * w for w in win)

    half = L // 2
    accum = [0.0] * (half + 1)
    n_seg = 0
    start = 0
    while start + L <= n:
        seg = x[start:start + L]
        m = sum(seg) / L
        windowed = [(seg[i] - m) * win[i] for i in range(L)]
        coeffs = dft(windowed)
        for k in range(half + 1):
            accum[k] += abs(coeffs[k]) ** 2 / win_power
        n_seg += 1
        start += step
    if n_seg == 0:
        raise ValueError("no full segments; reduce segment_length")
    freqs = [k / L for k in range(half + 1)]
    power = [accum[k] / n_seg for k in range(half + 1)]
    return freqs, power
