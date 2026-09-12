"""Robust scale and location estimators.

The standard deviation and the mean both have a breakdown point of zero: a single
extreme observation moves them arbitrarily far. These estimators degrade
gracefully instead.

  * ``median_absolute_deviation`` -- the median of ``|x - median(x)|``, scaled by
    1.4826 to be a consistent estimate of the standard deviation under normality.
    Breakdown point 50%.
  * ``interquartile_range`` -- Q3 - Q1, and its normal-consistent scale
    (divide by 1.349).
  * ``winsorize`` -- clip the tails to given percentiles (limits influence without
    dropping points).
  * ``trimmed_mean`` -- mean after discarding a fraction from each tail.

Pure standard library.
"""

_NORMAL_MAD_SCALE = 1.4826       # 1 / Phi^{-1}(0.75)
_NORMAL_IQR_SCALE = 1.349        # Phi^{-1}(0.75) - Phi^{-1}(0.25)


def _median(values):
    s = sorted(values)
    n = len(s)
    if n == 0:
        raise ValueError("empty sequence")
    mid = n // 2
    if n % 2 == 1:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def _quantile(sorted_values, q):
    """Linear-interpolation quantile of an already-sorted list (numpy default)."""
    n = len(sorted_values)
    if n == 1:
        return sorted_values[0]
    pos = q * (n - 1)
    lo = int(pos)
    frac = pos - lo
    if lo + 1 < n:
        return sorted_values[lo] * (1.0 - frac) + sorted_values[lo + 1] * frac
    return sorted_values[lo]


def median_absolute_deviation(x, scale=True):
    """Median absolute deviation of ``x``.

    ``MAD = median(|x_i - median(x)|)``. With ``scale=True`` (default) it is
    multiplied by 1.4826 so it consistently estimates the standard deviation for
    Gaussian data. Breakdown point 50%: up to half the data can be corrupted
    before it blows up. Zero for constant data.
    """
    if len(x) == 0:
        raise ValueError("need at least 1 point")
    med = _median(x)
    mad = _median([abs(v - med) for v in x])
    return mad * _NORMAL_MAD_SCALE if scale else mad


def interquartile_range(x, scale=False):
    """Interquartile range ``Q3 - Q1`` of ``x``.

    With ``scale=True`` it is divided by 1.349 to give a normal-consistent scale
    estimate (the IQR of a standard normal). Robust to outliers in the outer
    quartiles.
    """
    if len(x) < 2:
        raise ValueError("need at least 2 points")
    s = sorted(x)
    iqr = _quantile(s, 0.75) - _quantile(s, 0.25)
    return iqr / _NORMAL_IQR_SCALE if scale else iqr


def winsorize(x, limit=0.05):
    """Clip the tails of ``x`` to the ``limit`` / ``1 - limit`` quantiles.

    Returns a new list with values below the lower quantile raised to it and
    values above the upper quantile lowered to it -- bounding the influence of
    extremes without discarding observations. ``limit`` must be in ``[0, 0.5)``.
    """
    if not (0.0 <= limit < 0.5):
        raise ValueError("limit must be in [0, 0.5)")
    if len(x) == 0:
        raise ValueError("need at least 1 point")
    s = sorted(x)
    lo = _quantile(s, limit)
    hi = _quantile(s, 1.0 - limit)
    return [min(max(v, lo), hi) for v in x]


def trimmed_mean(x, proportion=0.1):
    """Mean of ``x`` after discarding a ``proportion`` fraction from each tail.

    ``proportion`` must be in ``[0, 0.5)``. With ``proportion = 0`` this is the
    ordinary mean; larger values give a more robust central estimate.
    """
    if not (0.0 <= proportion < 0.5):
        raise ValueError("proportion must be in [0, 0.5)")
    n = len(x)
    if n == 0:
        raise ValueError("need at least 1 point")
    s = sorted(x)
    k = int(n * proportion)
    kept = s[k:n - k] if n - k > k else s
    return sum(kept) / len(kept)
