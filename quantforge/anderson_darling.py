"""Anderson-Darling test for normality.

The Anderson-Darling statistic

    A^2 = -n - (1/n) sum_{i=1}^n (2i-1) [ln F(z_i) + ln(1 - F(z_{n+1-i}))]

weights the tails of the empirical-vs-normal CDF gap more heavily than
Kolmogorov-Smirnov, so it is the more powerful test for departures in the tails
(the ones that matter for financial returns). With the mean and standard deviation
estimated from the sample, ``A^2`` is adjusted for sample size and mapped to a
p-value by the D'Agostino-Stephens approximation. Pure standard library on top of
the normal CDF.
"""

import math

from .mathfns import norm_cdf


def anderson_darling_normal(values):
    """Anderson-Darling test that ``values`` are normal (mean/sd estimated).

    Returns ``(a2_star, p_value)`` where ``a2_star`` is the sample-size-adjusted
    statistic. A small p-value rejects normality; the test is especially sensitive to
    heavy tails and skew. Requires at least 8 observations for a meaningful p-value.
    """
    n = len(values)
    if n < 8:
        raise ValueError("need at least 8 observations")
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    if var <= 0.0:
        raise ValueError("zero-variance sample")
    sd = math.sqrt(var)
    z = sorted((v - mean) / sd for v in values)

    s = 0.0
    for i in range(n):
        fi = norm_cdf(z[i])
        fj = norm_cdf(z[n - 1 - i])
        # Guard the logs against 0/1 saturation in the far tails.
        fi = min(max(fi, 1e-300), 1.0 - 1e-16)
        fj = min(max(fj, 1e-300), 1.0 - 1e-16)
        s += (2 * (i + 1) - 1) * (math.log(fi) + math.log(1.0 - fj))
    a2 = -n - s / n
    # Small-sample adjustment (Stephens) for estimated parameters.
    a2_star = a2 * (1.0 + 4.0 / n - 25.0 / (n * n))
    return a2_star, _ad_pvalue(a2_star)


def _ad_pvalue(a2_star):
    """D'Agostino-Stephens p-value approximation for the adjusted A^2."""
    if a2_star < 0.2:
        return 1.0 - math.exp(-13.436 + 101.14 * a2_star - 223.73 * a2_star ** 2)
    if a2_star < 0.34:
        return 1.0 - math.exp(-8.318 + 42.796 * a2_star - 59.938 * a2_star ** 2)
    if a2_star < 0.6:
        return math.exp(0.9177 - 4.279 * a2_star - 1.38 * a2_star ** 2)
    if a2_star < 10.0:
        return math.exp(1.2937 - 5.709 * a2_star + 0.0186 * a2_star ** 2)
    return 0.0
