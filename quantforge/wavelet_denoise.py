"""Wavelet denoising by threshold shrinkage of Haar detail coefficients.

Donoho and Johnstone's wavelet shrinkage: transform a noisy signal to the wavelet
domain, shrink the small detail coefficients toward zero (they are mostly noise),
and invert. Smooth structure lives in a few large coefficients that survive; noise
is spread thinly across many small ones that get zeroed. The universal (VisuShrink)
threshold is ``sigma * sqrt(2 log n)``, with the noise scale ``sigma`` estimated
robustly from the finest detail level via the median absolute deviation. Pure
standard library; builds on :mod:`quantforge.wavelet`.
"""

import math

from .wavelet import haar_dwt, haar_idwt


def soft_threshold(x, lam):
    """Soft-threshold (shrink toward zero): ``sign(x) * max(|x| - lam, 0)``."""
    if x > lam:
        return x - lam
    if x < -lam:
        return x + lam
    return 0.0


def hard_threshold(x, lam):
    """Hard-threshold (keep or kill): ``x`` if ``|x| > lam`` else ``0``."""
    return x if abs(x) > lam else 0.0


def _median(values):
    s = sorted(values)
    n = len(s)
    if n == 0:
        return 0.0
    mid = n // 2
    if n % 2:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def mad_sigma(detail):
    """Robust noise-scale estimate from detail coefficients.

    ``sigma = median(|d|) / 0.6745`` -- the median absolute deviation rescaled to
    match the standard deviation of a Gaussian. Robust to the few large
    (signal-bearing) coefficients that would inflate a plain standard deviation.
    """
    if not detail:
        return 0.0
    return _median([abs(d) for d in detail]) / 0.6744897501960817


def universal_threshold(n, sigma):
    """VisuShrink universal threshold ``sigma * sqrt(2 log n)``."""
    if n < 2:
        return 0.0
    return sigma * math.sqrt(2.0 * math.log(n))


def wavelet_denoise(x, levels=None, mode="soft", threshold=None):
    """Denoise a signal by Haar wavelet shrinkage.

    Transforms ``x`` with :func:`haar_dwt`, shrinks every detail coefficient with the
    ``soft`` (default) or ``hard`` rule at the given ``threshold``, and inverts. If
    ``threshold`` is None the VisuShrink universal threshold is used, with the noise
    scale estimated by :func:`mad_sigma` from the finest detail level. The coarse
    approximation is left untouched (it carries the trend, not noise). Length must be
    a power of two. Returns the reconstructed, denoised signal.
    """
    if mode not in ("soft", "hard"):
        raise ValueError("mode must be 'soft' or 'hard'")
    approx, details = haar_dwt(x, levels)
    if threshold is None:
        sigma = mad_sigma(details[0])
        lam = universal_threshold(len(x), sigma)
    else:
        lam = threshold
    shrink = soft_threshold if mode == "soft" else hard_threshold
    new_details = [[shrink(c, lam) for c in level] for level in details]
    return haar_idwt(approx, new_details)
