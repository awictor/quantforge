"""Bland-Altman agreement analysis and Lin's concordance correlation coefficient.

Regression asks how two methods *relate*; agreement analysis asks how closely they
*match*. Two standard tools:

  * ``bland_altman`` -- the bias (mean difference) and the 95% limits of agreement
    ``bias +/- 1.96 * sd(diff)``, within which most method-to-method differences fall.
    It also returns the per-point means and differences for the classic Bland-Altman
    plot. A correlation can be high while agreement is poor (a constant offset), which
    this exposes directly.
  * ``concordance_correlation`` -- Lin's CCC, which combines precision (Pearson r) and
    accuracy (how close the best-fit line is to the 45-degree line of identity) into a
    single ``[-1, 1]`` agreement index. CCC = 1 only when every point lies on ``y = x``.

Pure standard library.
"""

import math


def bland_altman(x, y, k=1.96):
    """Bland-Altman agreement statistics for paired measurements.

    Returns a dict with the ``bias`` (mean of ``x - y``), ``sd`` of the differences,
    the ``lower`` and ``upper`` limits of agreement (``bias +/- k * sd``), and the
    per-point ``means`` and ``diffs`` for plotting. ``k`` defaults to 1.96 (95% limits).
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 paired points")
    diffs = [x[i] - y[i] for i in range(n)]
    means = [(x[i] + y[i]) / 2.0 for i in range(n)]
    bias = sum(diffs) / n
    var = sum((d - bias) ** 2 for d in diffs) / (n - 1)
    sd = math.sqrt(var)
    return {
        "bias": bias,
        "sd": sd,
        "lower": bias - k * sd,
        "upper": bias + k * sd,
        "means": means,
        "diffs": diffs,
    }


def concordance_correlation(x, y):
    """Lin's concordance correlation coefficient (CCC) of paired measurements.

    ``CCC = 2 s_xy / (s_x^2 + s_y^2 + (mx - my)^2)`` -- Pearson correlation penalized
    for any departure from the line of identity ``y = x``. Lies in ``[-1, 1]``: 1 only
    when the points fall exactly on ``y = x``, and it drops below the Pearson value
    whenever there is a systematic offset or scale difference.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have equal length")
    if n < 2:
        raise ValueError("need at least 2 paired points")
    mx = sum(x) / n
    my = sum(y) / n
    sx2 = sum((xi - mx) ** 2 for xi in x) / n
    sy2 = sum((yi - my) ** 2 for yi in y) / n
    sxy = sum((x[i] - mx) * (y[i] - my) for i in range(n)) / n
    denom = sx2 + sy2 + (mx - my) ** 2
    if denom == 0:
        return 1.0
    return 2.0 * sxy / denom
