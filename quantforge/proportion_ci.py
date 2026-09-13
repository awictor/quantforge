"""Confidence intervals for a binomial proportion.

Given ``k`` successes in ``n`` trials, four standard intervals for the underlying
success probability ``p``:

- ``wald`` -- the textbook normal approximation ``phat +/- z sqrt(phat(1-phat)/n)``
  (poor near 0 or 1 and for small ``n``),
- ``wilson`` -- the score interval, which stays inside ``[0, 1]`` and has far better
  small-sample coverage,
- ``agresti_coull`` -- add ``z^2/2`` pseudo-successes and failures, then a Wald
  interval on the adjusted counts, and
- ``clopper_pearson`` -- the exact interval by inverting the binomial CDF (never
  under-covers).

Each returns ``(low, high)``. Pure standard library on top of the binomial CDF and
the inverse error function.
"""

import math

from .special import erfinv
from .distributions import binomial_cdf


def _z(confidence):
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    # Two-sided normal critical value.
    return math.sqrt(2.0) * erfinv(confidence)


def _check(k, n):
    if n < 1:
        raise ValueError("n must be at least 1")
    if not (0 <= k <= n):
        raise ValueError("require 0 <= k <= n")


def wald_interval(k, n, confidence=0.95):
    """Normal-approximation (Wald) interval, clamped to ``[0, 1]``."""
    _check(k, n)
    z = _z(confidence)
    phat = k / n
    half = z * math.sqrt(phat * (1.0 - phat) / n)
    return max(0.0, phat - half), min(1.0, phat + half)


def wilson_interval(k, n, confidence=0.95):
    """Wilson score interval -- stays in ``[0, 1]`` with good small-sample coverage."""
    _check(k, n)
    z = _z(confidence)
    phat = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (phat + z2 / (2.0 * n)) / denom
    half = (z * math.sqrt(phat * (1.0 - phat) / n + z2 / (4.0 * n * n))) / denom
    return max(0.0, center - half), min(1.0, center + half)


def agresti_coull_interval(k, n, confidence=0.95):
    """Agresti-Coull interval: a Wald interval on ``z^2``-adjusted counts."""
    _check(k, n)
    z = _z(confidence)
    z2 = z * z
    n_adj = n + z2
    p_adj = (k + z2 / 2.0) / n_adj
    half = z * math.sqrt(p_adj * (1.0 - p_adj) / n_adj)
    return max(0.0, p_adj - half), min(1.0, p_adj + half)


def clopper_pearson_interval(k, n, confidence=0.95):
    """Exact Clopper-Pearson interval by inverting the binomial CDF.

    The lower limit is the ``p`` with ``P(X >= k) = alpha/2`` and the upper limit the
    ``p`` with ``P(X <= k) = alpha/2`` (``alpha = 1 - confidence``); the boundary
    cases ``k = 0`` and ``k = n`` give a one-sided interval. Guaranteed to cover at
    least ``confidence`` of the time -- conservative but never under-covering.
    """
    _check(k, n)
    alpha = 1.0 - confidence
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")

    # Lower: solve P(X >= k | p) = alpha/2, i.e. 1 - binomial_cdf(k-1, n, p) = alpha/2.
    if k == 0:
        low = 0.0
    else:
        low = _bisect(lambda p: 1.0 - binomial_cdf(k - 1, n, p) - alpha / 2.0, 0.0, 1.0)
    # Upper: solve P(X <= k | p) = alpha/2, i.e. binomial_cdf(k, n, p) = alpha/2.
    if k == n:
        high = 1.0
    else:
        high = _bisect(lambda p: binomial_cdf(k, n, p) - alpha / 2.0, 0.0, 1.0)
    return low, high


def _bisect(f, lo, hi):
    """Bisection root of a monotone ``f`` on ``[lo, hi]`` (200 iterations)."""
    flo = f(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if (fmid < 0.0) == (flo < 0.0):
            lo, flo = mid, fmid
        else:
            hi = mid
        if hi - lo < 1e-14:
            break
    return 0.5 * (lo + hi)
