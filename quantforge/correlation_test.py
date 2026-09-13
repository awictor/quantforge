"""Significance test and confidence interval for the Pearson correlation.

Given paired samples ``x`` and ``y``, the sample correlation ``r`` is tested for
departure from zero with the exact-under-normality t statistic

    t = r sqrt((n - 2) / (1 - r^2)),   df = n - 2,

and a confidence interval is built with Fisher's variance-stabilizing z-transform
``z = atanh(r)``, which is approximately normal with standard error
``1 / sqrt(n - 3)``. Pure standard library on top of the t distribution and the
inverse error function.
"""

import math

from .student_t import t_cdf
from .special import erfinv


def _mean(v):
    return sum(v) / len(v)


def pearson_r(x, y):
    """Pearson product-moment correlation coefficient of paired samples.

    ``r = cov(x, y) / (sd(x) sd(y))`` in ``[-1, 1]``. Raises if either sample has
    zero variance.
    """
    n = len(x)
    if n < 2 or len(y) != n:
        raise ValueError("x and y must be equal-length with at least two points")
    mx, my = _mean(x), _mean(y)
    sxy = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    sxx = sum((x[i] - mx) ** 2 for i in range(n))
    syy = sum((y[i] - my) ** 2 for i in range(n))
    if sxx <= 0.0 or syy <= 0.0:
        raise ValueError("each sample must have non-zero variance")
    return sxy / math.sqrt(sxx * syy)


def pearson_correlation_test(x, y, confidence=0.95):
    """Test ``H0: rho = 0`` and give a Fisher-z confidence interval for ``rho``.

    Returns a dict with ``r`` (the sample correlation), ``t_stat`` and ``df`` of the
    two-sided t-test, ``p_value``, and ``conf_int`` ``[low, high]`` at ``confidence``
    from Fisher's z-transform. The p-value is small when the correlation is unlikely
    to be zero; the interval is clamped to ``[-1, 1]`` and, at ``|r| = 1``, collapses
    to the point.
    """
    n = len(x)
    if n < 3:
        raise ValueError("need at least three points for a confidence interval")
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    r = pearson_r(x, y)
    df = n - 2

    if abs(r) >= 1.0:
        return {"r": r, "t_stat": math.inf * (1.0 if r > 0 else -1.0), "df": df,
                "p_value": 0.0, "conf_int": [r, r]}

    t_stat = r * math.sqrt(df / (1.0 - r * r))
    p_value = 2.0 * (1.0 - t_cdf(abs(t_stat), df))

    # Fisher z-transform interval.
    z = math.atanh(r)
    se = 1.0 / math.sqrt(n - 3)
    # Normal quantile at the two-sided level via erfinv.
    z_crit = math.sqrt(2.0) * erfinv(confidence)
    lo = math.tanh(z - z_crit * se)
    hi = math.tanh(z + z_crit * se)
    return {"r": r, "t_stat": t_stat, "df": df, "p_value": p_value,
            "conf_int": [max(-1.0, lo), min(1.0, hi)]}
