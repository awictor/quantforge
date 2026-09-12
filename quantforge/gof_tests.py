"""Goodness-of-fit tests: Jarque-Bera (with p-value) and two-sample KS.

  * ``jarque_bera_test`` -- the Jarque-Bera normality statistic together with its
    chi-square(2) p-value (the bare statistic lives in ``perfmetrics``). Large
    statistic / small p rejects normality.
  * ``ks_two_sample`` -- the two-sample Kolmogorov-Smirnov statistic ``D``, the
    largest gap between the two empirical CDFs, with the asymptotic p-value. It
    detects any distributional difference (location, scale, shape), not just a
    difference in means.

Pure standard library.
"""

import math

from .perfmetrics import jarque_bera
from .serial_correlation import _chi2_sf


def jarque_bera_test(returns):
    """Jarque-Bera normality test: ``(statistic, p_value)``.

    The statistic is chi-square(2) under the normal null, so the p-value is its
    upper-tail probability. A normal sample gives a small statistic and a large
    p-value; a heavy-tailed or skewed sample gives a large statistic and a small
    p-value.
    """
    jb = jarque_bera(returns)
    return jb, _chi2_sf(jb, 2)


def _ecdf_at(sorted_sample, x):
    """Empirical CDF of ``sorted_sample`` evaluated at ``x`` (fraction <= x)."""
    lo, hi = 0, len(sorted_sample)
    while lo < hi:
        mid = (lo + hi) // 2
        if sorted_sample[mid] <= x:
            lo = mid + 1
        else:
            hi = mid
    return lo / len(sorted_sample)


def _ks_pvalue(d, en):
    """Asymptotic Kolmogorov distribution p-value for statistic ``d``.

    ``P(D > d) = 2 sum_{k=1}^inf (-1)^{k-1} exp(-2 k^2 (en d)^2)`` where ``en`` is
    the effective sample size. Clamped to [0, 1].
    """
    t = (en + 0.12 + 0.11 / en) * d
    if t <= 0.0:
        return 1.0
    s = 0.0
    for k in range(1, 101):
        term = 2.0 * (-1.0) ** (k - 1) * math.exp(-2.0 * k * k * t * t)
        s += term
        if abs(term) < 1e-12:
            break
    return max(0.0, min(1.0, s))


def ks_two_sample(a, b):
    """Two-sample Kolmogorov-Smirnov test: ``(D, p_value)``.

    ``D`` is the maximum absolute difference between the two empirical CDFs,
    evaluated at every observed point. The p-value uses the asymptotic Kolmogorov
    distribution with the effective sample size ``sqrt(n m / (n + m))``. A small p
    rejects "the two samples come from the same distribution".
    """
    n, m = len(a), len(b)
    if n < 1 or m < 1:
        raise ValueError("both samples must be non-empty")
    sa = sorted(a)
    sb = sorted(b)
    d = 0.0
    for x in sa + sb:
        gap = abs(_ecdf_at(sa, x) - _ecdf_at(sb, x))
        if gap > d:
            d = gap
    en = math.sqrt(n * m / (n + m))
    return d, _ks_pvalue(d, en)
