"""D'Agostino-Pearson K^2 omnibus test of normality.

Combines two shape diagnostics into one normality test: the sample *skewness* (should
be 0 for a normal) and *kurtosis* (should be 3). Each is transformed to an approximately
standard-normal Z score -- D'Agostino's skewness transform and Anscombe-Glynn's kurtosis
transform -- and the sum of their squares is the K^2 statistic, chi-square with 2
degrees of freedom under normality. It catches both skew and heavy/light tails, unlike a
single-moment test. Needs a moderate sample (>= 20). Pure standard library.
"""

import math

from .distributions import chi2_cdf


def dagostino_k2(values):
    """D'Agostino-Pearson K^2 omnibus normality test.

    Returns a dict with the ``k2`` statistic (chi-square, 2 df), its ``p_value``, and
    the component ``z_skew`` and ``z_kurt`` standard scores. A small p-value rejects
    normality; the components show whether skew, tails, or both drive it. Needs at least
    20 observations for the transforms to be reliable.
    """
    n = len(values)
    if n < 8:
        raise ValueError("need at least 8 observations (>= 20 recommended)")
    m = sum(values) / n
    m2 = sum((x - m) ** 2 for x in values) / n
    if m2 <= 0:
        raise ValueError("zero variance; test undefined")
    m3 = sum((x - m) ** 3 for x in values) / n
    m4 = sum((x - m) ** 4 for x in values) / n
    b1 = m3 / m2 ** 1.5          # sample skewness
    b2 = m4 / (m2 * m2)          # sample kurtosis (normal -> 3)

    # D'Agostino (1970) skewness transform to a standard normal Z.
    y = b1 * math.sqrt((n + 1) * (n + 3) / (6.0 * (n - 2)))
    beta2 = (3.0 * (n * n + 27 * n - 70) * (n + 1) * (n + 3)
             / ((n - 2) * (n + 5) * (n + 7) * (n + 9)))
    w2 = -1.0 + math.sqrt(2.0 * (beta2 - 1.0))
    delta = 1.0 / math.sqrt(0.5 * math.log(w2))
    alpha = math.sqrt(2.0 / (w2 - 1.0))
    yy = y / alpha if alpha != 0 else y
    z_skew = delta * math.log(yy + math.sqrt(yy * yy + 1.0))

    # Anscombe-Glynn kurtosis transform.
    mean_b2 = 3.0 * (n - 1) / (n + 1)
    var_b2 = 24.0 * n * (n - 2) * (n - 3) / ((n + 1) ** 2 * (n + 3) * (n + 5))
    x = (b2 - mean_b2) / math.sqrt(var_b2)
    sqrt_beta1 = (6.0 * (n * n - 5 * n + 2) / ((n + 7) * (n + 9))
                  * math.sqrt(6.0 * (n + 3) * (n + 5) / (n * (n - 2) * (n - 3))))
    a = 6.0 + 8.0 / sqrt_beta1 * (2.0 / sqrt_beta1
                                  + math.sqrt(1.0 + 4.0 / (sqrt_beta1 ** 2)))
    term = (1.0 - 2.0 / a) / (1.0 + x * math.sqrt(2.0 / (a - 4.0)))
    z_kurt = ((1.0 - 2.0 / (9.0 * a)) - term ** (1.0 / 3.0)) / math.sqrt(2.0 / (9.0 * a))

    k2 = z_skew * z_skew + z_kurt * z_kurt
    p = 1.0 - chi2_cdf(k2, 2)
    return {"k2": k2, "p_value": p, "z_skew": z_skew, "z_kurt": z_kurt}
