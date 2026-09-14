"""Robust scale estimators: Qn, Sn, and the biweight midvariance.

High-breakdown alternatives to the standard deviation for spread. The Rousseeuw-Croux Qn
and Sn estimators tolerate up to 50% contamination and, unlike the MAD, are not tied to
symmetry (Sn) or need a location estimate. The biweight midvariance smoothly downweights
outliers. Each carries the consistency constant that makes it an unbiased estimate of the
normal standard deviation on clean Gaussian data. Pure standard library.
"""

import math


def qn_scale(x):
    """Rousseeuw-Croux Qn scale estimator (50% breakdown, location-free).

    The (roughly) first quartile of the pairwise distances ``|x_i - x_j|`` (``i < j``),
    scaled by ``2.2219`` for asymptotic consistency with the normal standard deviation.
    Needs at least two points; robust to up to half the data being outliers. No
    finite-sample correction is applied, so small-sample values differ slightly from
    implementations that include one.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least two values")
    diffs = []
    for i in range(n):
        for j in range(i + 1, n):
            diffs.append(abs(x[i] - x[j]))
    diffs.sort()
    h = n // 2 + 1
    k = h * (h - 1) // 2                      # target order statistic (1-indexed)
    q = diffs[k - 1]
    return 2.2219 * q


def sn_scale(x):
    """Rousseeuw-Croux Sn scale estimator (50% breakdown, no location needed).

    ``1.1926 * median_i( median_j |x_i - x_j| )`` -- the outer/inner medians give a robust
    spread that, unlike the MAD, does not assume a symmetric distribution. Needs at least
    two points.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least two values")
    med_of = []
    for i in range(n):
        dists = sorted(abs(x[i] - x[j]) for j in range(n))
        # "high median": order statistic at index floor(n/2) (Rousseeuw-Croux convention).
        med_of.append(dists[n // 2])
    med_of.sort()
    # Low median of the per-point high medians.
    c = med_of[(n - 1) // 2]
    return 1.1926 * c


def biweight_midvariance(x, c=9.0):
    """Biweight midvariance: a robust variance that smoothly downweights outliers.

    Returns the square root (a robust standard deviation). Points more than ``c`` MADs
    from the median get zero weight; the tuning constant ``c=9`` gives ~87% efficiency at
    the normal. Reduces to a near-standard-deviation on clean Gaussian data.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least two values")
    s = sorted(x)
    med = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
    devs = sorted(abs(v - med) for v in x)
    mad = devs[n // 2] if n % 2 else 0.5 * (devs[n // 2 - 1] + devs[n // 2])
    if mad == 0.0:
        return 0.0
    num = 0.0
    den = 0.0
    for v in x:
        u = (v - med) / (c * mad)
        if abs(u) < 1.0:
            w = (1.0 - u * u)
            num += (v - med) ** 2 * w ** 4
            den += w * (1.0 - 5.0 * u * u)
    if den == 0.0:
        return 0.0
    return math.sqrt(n * num) / abs(den)
