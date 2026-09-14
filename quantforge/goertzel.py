"""Goertzel algorithm: efficient single-frequency DFT.

When you need the strength of one (or a few) known frequencies rather than the whole
spectrum, the Goertzel algorithm computes a single DFT bin in ``O(n)`` with one real
recurrence -- far cheaper than an FFT and the standard method for tone detection (DTMF,
pilot tones, a known harmonic). Returns the complex DFT coefficient at that bin and its
power. Pure standard library.
"""

import math


def goertzel(x, k):
    """DFT coefficient ``X[k]`` of ``x`` at integer bin ``k`` by the Goertzel recurrence.

    Returns the complex ``X[k] = sum_n x[n] exp(-2 pi i k n / N)`` -- identical to the
    ``k``-th FFT output, computed with a single ``O(n)`` real recurrence. Integer ``k``
    in ``[0, N)``.
    """
    n = len(x)
    if n == 0:
        raise ValueError("need at least one sample")
    if not (0 <= k < n):
        raise ValueError("bin k must be in [0, N)")
    w = 2.0 * math.pi * k / n
    coeff = 2.0 * math.cos(w)
    s1 = 0.0
    s2 = 0.0
    for sample in x:
        s0 = sample + coeff * s1 - s2
        s2 = s1
        s1 = s0
    # X[k] = e^{jw} (s1 - e^{-jw} s2); the leading e^{jw} sets the correct phase so
    # this matches the k-th DFT/FFT coefficient exactly, not just in magnitude.
    return complex(math.cos(w), math.sin(w)) * (s1 - complex(math.cos(w), -math.sin(w)) * s2)


def goertzel_power(x, k):
    """Power ``|X[k]|^2`` at bin ``k`` (the Goertzel magnitude-squared).

    The efficient tone-detection quantity: large when a frequency near bin ``k`` is
    present, small otherwise. Avoids the final trig of :func:`goertzel`.
    """
    n = len(x)
    if n == 0:
        raise ValueError("need at least one sample")
    if not (0 <= k < n):
        raise ValueError("bin k must be in [0, N)")
    w = 2.0 * math.pi * k / n
    coeff = 2.0 * math.cos(w)
    s1 = 0.0
    s2 = 0.0
    for sample in x:
        s0 = sample + coeff * s1 - s2
        s2 = s1
        s1 = s0
    return s1 * s1 + s2 * s2 - coeff * s1 * s2
