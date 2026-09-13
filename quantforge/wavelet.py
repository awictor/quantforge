"""Haar discrete wavelet transform (multilevel) and energy decomposition.

The Haar wavelet is the simplest orthonormal wavelet: at each level it replaces
adjacent pairs by their scaled sum (approximation) and difference (detail),

    a = (x_{2i} + x_{2i+1}) / sqrt(2),   d = (x_{2i} - x_{2i+1}) / sqrt(2),

then recurses on the approximation. This gives a multiresolution view -- coarse
trend plus detail coefficients at each scale -- and, being orthonormal, preserves
total energy (Parseval) and reconstructs exactly. Length must be a power of two.
Pure standard library.
"""

import math

_SQRT2 = math.sqrt(2.0)


def haar_dwt(x, levels=None):
    """Multilevel Haar wavelet transform.

    Returns ``(approx, details)`` where ``approx`` is the final coarse-approximation
    list and ``details`` is a list (finest level first) of the detail-coefficient
    lists at each level. ``levels`` defaults to the maximum ``log2(len(x))``. Length
    must be a power of two.
    """
    n = len(x)
    if n == 0 or (n & (n - 1)) != 0:
        raise ValueError("length must be a positive power of two")
    max_levels = n.bit_length() - 1
    if levels is None:
        levels = max_levels
    if not (1 <= levels <= max_levels):
        raise ValueError("levels must be in [1, log2(n)]")
    a = list(x)
    details = []
    for _ in range(levels):
        m = len(a) // 2
        approx = [(a[2 * i] + a[2 * i + 1]) / _SQRT2 for i in range(m)]
        detail = [(a[2 * i] - a[2 * i + 1]) / _SQRT2 for i in range(m)]
        details.append(detail)
        a = approx
    return a, details


def haar_idwt(approx, details):
    """Invert :func:`haar_dwt`, reconstructing the original signal.

    Takes the coarse approximation and the per-level detail lists (finest first) and
    returns the reconstructed series. Exact up to floating error.
    """
    a = list(approx)
    for detail in reversed(details):
        m = len(a)
        out = [0.0] * (2 * m)
        for i in range(m):
            out[2 * i] = (a[i] + detail[i]) / _SQRT2
            out[2 * i + 1] = (a[i] - detail[i]) / _SQRT2
        a = out
    return a


def wavelet_energy(x, levels=None):
    """Fraction of signal energy in each Haar detail level and the coarse approximation.

    Returns a dict with ``detail`` (a list of energy fractions, finest level first)
    and ``approx`` (the coarse-approximation energy fraction). The fractions sum to
    one because the Haar transform is orthonormal (Parseval). A smooth series
    concentrates energy in the approximation; a noisy one spreads it into the fine
    details.
    """
    total = sum(v * v for v in x)
    if total <= 0.0:
        raise ValueError("signal has zero energy")
    a, details = haar_dwt(x, levels)
    detail_frac = [sum(d * d for d in lvl) / total for lvl in details]
    approx_frac = sum(v * v for v in a) / total
    return {"detail": detail_frac, "approx": approx_frac}
