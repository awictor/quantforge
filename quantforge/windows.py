"""Window functions for spectral analysis.

Cutting a finite segment out of a signal is an implicit multiplication by a rectangular
window, whose sharp edges leak energy across the spectrum (sidelobes). Tapering the
segment to zero at the ends with a smooth window trades a little main-lobe width for far
lower sidelobes -- the standard pre-processing before an FFT/periodogram. Provides the
common windows (Hann, Hamming, Blackman, Bartlett, rectangular) and an ``apply_window``
helper. Pure standard library.
"""

import math


def hann(n):
    """Hann (raised-cosine) window of length ``n``: ``0.5(1 - cos(2 pi k/(n-1)))``."""
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1:
        return [1.0]
    return [0.5 * (1.0 - math.cos(2.0 * math.pi * k / (n - 1))) for k in range(n)]


def hamming(n):
    """Hamming window: ``0.54 - 0.46 cos(2 pi k/(n-1))`` (lower first sidelobe than Hann)."""
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1:
        return [1.0]
    return [0.54 - 0.46 * math.cos(2.0 * math.pi * k / (n - 1)) for k in range(n)]


def blackman(n):
    """Blackman window: three-term cosine, very low sidelobes at a wider main lobe."""
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1:
        return [1.0]
    a0, a1, a2 = 0.42, 0.5, 0.08
    return [a0 - a1 * math.cos(2.0 * math.pi * k / (n - 1))
            + a2 * math.cos(4.0 * math.pi * k / (n - 1)) for k in range(n)]


def bartlett(n):
    """Bartlett (triangular) window, zero at both ends."""
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1:
        return [1.0]
    m = (n - 1) / 2.0
    return [1.0 - abs((k - m) / m) for k in range(n)]


def rectangular(n):
    """Rectangular (boxcar) window: all ones -- no tapering."""
    if n < 1:
        raise ValueError("n must be >= 1")
    return [1.0] * n


def apply_window(x, window="hann"):
    """Multiply signal ``x`` by a named window (or a precomputed window list).

    ``window`` is ``"hann"``, ``"hamming"``, ``"blackman"``, ``"bartlett"``,
    ``"rectangular"``, or a list of the same length as ``x``. Returns the tapered
    signal.
    """
    n = len(x)
    if isinstance(window, str):
        table = {"hann": hann, "hamming": hamming, "blackman": blackman,
                 "bartlett": bartlett, "rectangular": rectangular}
        if window not in table:
            raise ValueError("unknown window name")
        w = table[window](n)
    else:
        w = list(window)
        if len(w) != n:
            raise ValueError("window length must match signal")
    return [x[i] * w[i] for i in range(n)]
