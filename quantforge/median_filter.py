"""Nonlinear order-statistic filters: median, rank, and Hampel.

Linear filters smear sharp features and are dragged by outliers. Order-statistic filters
sort the samples in a sliding window and pick one by rank, so they reject impulsive
spikes while preserving edges -- the classic tool for de-spiking sensor data and for
edge-preserving smoothing. ``median_filter`` takes the window median, ``rank_filter``
any percentile, and ``hampel_filter`` replaces only the points that sit too many robust
standard deviations from the local median (leaving clean data untouched). Pure standard
library.
"""

import math


def _window_bounds(i, n, half):
    lo = max(0, i - half)
    hi = min(n, i + half + 1)
    return lo, hi


def median_filter(x, window):
    """Sliding-window median of ``x`` with an odd ``window`` length.

    Replaces each point by the median of the ``window`` samples centered on it (the
    window is clipped at the ends). Removes impulsive spikes while keeping step edges
    sharp, unlike a moving average. Returns a list the same length as ``x``.
    """
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    if window < 1 or window % 2 == 0:
        raise ValueError("window must be a positive odd integer")
    half = window // 2
    out = []
    for i in range(n):
        lo, hi = _window_bounds(i, n, half)
        seg = sorted(x[lo:hi])
        m = len(seg)
        out.append(seg[m // 2] if m % 2 else 0.5 * (seg[m // 2 - 1] + seg[m // 2]))
    return out


def rank_filter(x, window, percentile):
    """Sliding-window order-statistic filter at a given ``percentile`` (0-100).

    ``percentile=50`` is the median; ``0`` a min filter, ``100`` a max filter. Picks the
    order statistic nearest that percentile in each centered window. Returns a list the
    same length as ``x``.
    """
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    if window < 1 or window % 2 == 0:
        raise ValueError("window must be a positive odd integer")
    if not (0.0 <= percentile <= 100.0):
        raise ValueError("percentile must be in [0, 100]")
    half = window // 2
    out = []
    for i in range(n):
        lo, hi = _window_bounds(i, n, half)
        seg = sorted(x[lo:hi])
        idx = int(round((percentile / 100.0) * (len(seg) - 1)))
        out.append(seg[idx])
    return out


def hampel_filter(x, window=7, n_sigmas=3.0):
    """Hampel outlier filter: replace points far from the local median with that median.

    In each centered window of length ``2*window+1`` it computes the median and the
    median absolute deviation (MAD), scales the MAD to a robust standard deviation
    (``1.4826 * MAD``), and replaces the center point only if it lies more than
    ``n_sigmas`` robust deviations away. Clean data passes through untouched. Returns
    ``(filtered, outlier_indices)``.
    """
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    if window < 1:
        raise ValueError("window must be >= 1")
    if n_sigmas <= 0.0:
        raise ValueError("n_sigmas must be positive")
    k = 1.4826            # MAD -> sigma for normally-distributed data
    out = list(x)
    outliers = []
    for i in range(n):
        lo, hi = _window_bounds(i, n, window)
        seg = sorted(x[lo:hi])
        m = len(seg)
        med = seg[m // 2] if m % 2 else 0.5 * (seg[m // 2 - 1] + seg[m // 2])
        devs = sorted(abs(v - med) for v in x[lo:hi])
        mad = devs[m // 2] if m % 2 else 0.5 * (devs[m // 2 - 1] + devs[m // 2])
        sigma = k * mad
        if sigma > 0.0 and abs(x[i] - med) > n_sigmas * sigma:
            out[i] = med
            outliers.append(i)
    return out, outliers
