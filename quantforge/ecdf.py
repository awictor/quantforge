"""Empirical CDF, sample quantiles, and Q-Q comparison.

The building blocks under the distribution tests: the empirical CDF ``F_n(x)`` (the
fraction of sample points at or below ``x``), the sample quantile (inverse CDF, with the
common interpolation conventions), and the quantile-quantile pairing used to compare a
sample against another sample or a theoretical distribution. Pure standard library.
"""


def ecdf(data, x):
    """Empirical CDF ``F_n(x) = #{data_i <= x} / n`` at point(s) ``x``.

    ``x`` may be a scalar (returns a float in ``[0, 1]``) or an iterable (returns a
    list). Right-continuous step function.
    """
    n = len(data)
    if n == 0:
        raise ValueError("need at least one point")
    s = sorted(data)

    def f(xi):
        # Count of values <= xi via binary search (rightmost position).
        lo, hi = 0, n
        while lo < hi:
            mid = (lo + hi) // 2
            if s[mid] <= xi:
                lo = mid + 1
            else:
                hi = mid
        return lo / n

    if hasattr(x, "__iter__"):
        return [f(xi) for xi in x]
    return f(x)


def quantile(data, p, method="linear"):
    """Sample quantile at probability ``p`` in ``[0, 1]``.

    ``method`` selects the interpolation on the order statistics: ``"linear"`` (the
    NumPy default, ``(n-1)p`` position), ``"lower"``, ``"higher"``, or ``"nearest"``.
    ``p`` may be a scalar or iterable. Matches the standard percentile conventions.
    """
    n = len(data)
    if n == 0:
        raise ValueError("need at least one point")
    s = sorted(data)

    def q(pp):
        if not (0.0 <= pp <= 1.0):
            raise ValueError("p must be in [0, 1]")
        idx = pp * (n - 1)
        lo = int(idx)
        frac = idx - lo
        if method == "linear":
            if lo + 1 < n:
                return s[lo] * (1 - frac) + s[lo + 1] * frac
            return s[lo]
        if method == "lower":
            return s[lo]
        if method == "higher":
            return s[min(lo + 1, n - 1)] if frac > 0 else s[lo]
        if method == "nearest":
            return s[lo + 1] if (frac > 0.5 and lo + 1 < n) else s[lo]
        raise ValueError("method must be linear/lower/higher/nearest")

    if hasattr(p, "__iter__"):
        return [q(pp) for pp in p]
    return q(p)


def qq_points(sample, reference):
    """Quantile-quantile pairing of two samples for a Q-Q plot.

    Returns a list of ``(reference_quantile, sample_quantile)`` pairs evaluated at the
    plotting positions ``(i - 0.5)/m`` of the smaller sample size ``m``. On the y = x
    line the two samples share a distribution; a slope or curvature departure reveals a
    scale or shape difference. ``reference`` may be another sample.
    """
    m = min(len(sample), len(reference))
    if m == 0:
        raise ValueError("both inputs must be non-empty")
    probs = [(i + 0.5) / m for i in range(m)]
    sq = quantile(sample, probs)
    rq = quantile(reference, probs)
    return list(zip(rq, sq))
