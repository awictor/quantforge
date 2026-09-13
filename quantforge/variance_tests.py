"""Tests for equal variance across k groups (homogeneity of variance).

Many procedures (ANOVA, pooled t) assume the groups share a variance; these test that:

  * ``levene_test`` -- one-way ANOVA on the *absolute deviations* from each group's
    center. With ``center="median"`` it is the Brown-Forsythe variant, robust to
    non-normal data; ``center="mean"`` is the original Levene. F-distributed.
  * ``bartlett_test`` -- a likelihood-ratio test, more powerful under normality but
    sensitive to departures from it. Chi-square distributed.

Both return the statistic and a p-value. Pure standard library.
"""

import math

from .distributions import f_cdf, chi2_cdf


def _mean(v):
    return sum(v) / len(v)


def _median(v):
    s = sorted(v)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else 0.5 * (s[mid - 1] + s[mid])


def levene_test(*groups, center="median"):
    """Levene / Brown-Forsythe test for equal variance across groups.

    ``center="median"`` (default) is the robust Brown-Forsythe form; ``"mean"`` is the
    original Levene. Returns a dict with the ``statistic`` (F), ``df`` ``(k-1, N-k)`` and
    the ``p_value``. A small p-value rejects equal variance.
    """
    gs = [list(g) for g in groups if len(g) > 0]
    k = len(gs)
    if k < 2:
        raise ValueError("need at least 2 non-empty groups")
    ctr = _median if center == "median" else _mean if center == "mean" else None
    if ctr is None:
        raise ValueError("center must be 'median' or 'mean'")

    # z_ij = |x_ij - center_i|.
    z = [[abs(x - ctr(g)) for x in g] for g in gs]
    ns = [len(g) for g in gs]
    N = sum(ns)
    zbar = [_mean(zi) for zi in z]
    zbar_all = sum(sum(zi) for zi in z) / N

    numer = (N - k) * sum(ns[i] * (zbar[i] - zbar_all) ** 2 for i in range(k))
    denom = (k - 1) * sum((z[i][j] - zbar[i]) ** 2 for i in range(k) for j in range(ns[i]))
    if denom <= 0:
        raise ValueError("zero within-group deviation; test undefined")
    W = numer / denom
    df1, df2 = k - 1, N - k
    p = 1.0 - f_cdf(W, df1, df2)
    return {"statistic": W, "df": (df1, df2), "p_value": p}


def bartlett_test(*groups):
    """Bartlett's test for equal variance (likelihood ratio, chi-square).

    Returns a dict with the ``statistic``, ``df`` (``k-1``) and the ``p_value``. More
    powerful than Levene under normality but sensitive to non-normal tails. Requires at
    least two observations per group.
    """
    gs = [list(g) for g in groups if len(g) > 1]
    k = len(gs)
    if k < 2:
        raise ValueError("need at least 2 groups with >= 2 observations each")
    ns = [len(g) for g in gs]
    N = sum(ns)
    # Sample variances (n-1).
    variances = []
    for g in gs:
        m = _mean(g)
        variances.append(sum((x - m) ** 2 for x in g) / (len(g) - 1))
    if any(v <= 0 for v in variances):
        raise ValueError("a group has zero variance; test undefined")
    sp2 = sum((ns[i] - 1) * variances[i] for i in range(k)) / (N - k)
    numer = (N - k) * math.log(sp2) - sum((ns[i] - 1) * math.log(variances[i])
                                          for i in range(k))
    c = 1.0 + (sum(1.0 / (ns[i] - 1) for i in range(k)) - 1.0 / (N - k)) / (3 * (k - 1))
    stat = numer / c
    df = k - 1
    p = 1.0 - chi2_cdf(stat, df)
    return {"statistic": stat, "df": df, "p_value": p}
