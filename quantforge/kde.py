"""Gaussian kernel density estimation.

A histogram's shape depends on bin edges; a kernel density estimate is a smooth,
edge-free density built by summing a Gaussian bump over each data point,

    f(x) = (1 / (n h)) sum_i phi((x - x_i) / h),

with bandwidth ``h`` chosen by a rule of thumb (Silverman or Scott) from the sample
size and spread. The estimate integrates to one and converges to the true density as
``n`` grows. Provides bandwidth selection, pointwise evaluation, and a callable factory.
Pure standard library.
"""

import math

_INV_SQRT_2PI = 1.0 / math.sqrt(2.0 * math.pi)


def _std(data):
    n = len(data)
    m = sum(data) / n
    var = sum((x - m) ** 2 for x in data) / (n - 1)
    return math.sqrt(var)


def _iqr(data):
    s = sorted(data)
    n = len(s)

    def q(p):
        idx = p * (n - 1)
        lo = int(idx)
        frac = idx - lo
        if lo + 1 < n:
            return s[lo] * (1 - frac) + s[lo + 1] * frac
        return s[lo]

    return q(0.75) - q(0.25)


def silverman_bandwidth(data):
    """Silverman's rule-of-thumb bandwidth.

    ``h = 0.9 * min(std, IQR/1.34) * n^{-1/5}`` -- robust to mild non-normality via the
    IQR term. The standard default for a unimodal, roughly-Gaussian sample.
    """
    n = len(data)
    if n < 2:
        raise ValueError("need at least 2 points")
    sd = _std(data)
    iqr = _iqr(data)
    spread = min(sd, iqr / 1.34) if iqr > 0 else sd
    if spread <= 0:
        spread = sd if sd > 0 else 1.0
    return 0.9 * spread * n ** (-0.2)


def scott_bandwidth(data):
    """Scott's rule-of-thumb bandwidth ``h = std * n^{-1/5}``."""
    n = len(data)
    if n < 2:
        raise ValueError("need at least 2 points")
    sd = _std(data)
    return (sd if sd > 0 else 1.0) * n ** (-0.2)


def kde(data, x, bandwidth=None, rule="silverman"):
    """Evaluate the Gaussian KDE of ``data`` at point(s) ``x``.

    ``bandwidth`` overrides the rule (``"silverman"`` or ``"scott"``). ``x`` may be a
    scalar (returns a float) or an iterable (returns a list). The estimate is
    non-negative everywhere and integrates to one.
    """
    n = len(data)
    if n < 2:
        raise ValueError("need at least 2 points")
    if bandwidth is None:
        bandwidth = (silverman_bandwidth(data) if rule == "silverman"
                     else scott_bandwidth(data) if rule == "scott"
                     else None)
        if bandwidth is None:
            raise ValueError("rule must be 'silverman' or 'scott'")
    if bandwidth <= 0:
        raise ValueError("bandwidth must be positive")

    def density(xi):
        s = 0.0
        for d in data:
            u = (xi - d) / bandwidth
            s += _INV_SQRT_2PI * math.exp(-0.5 * u * u)
        return s / (n * bandwidth)

    if hasattr(x, "__iter__"):
        return [density(xi) for xi in x]
    return density(x)


def kde_function(data, bandwidth=None, rule="silverman"):
    """Return a callable density estimator ``f(x)`` for ``data`` (bandwidth fixed once)."""
    n = len(data)
    if n < 2:
        raise ValueError("need at least 2 points")
    if bandwidth is None:
        bandwidth = (silverman_bandwidth(data) if rule == "silverman"
                     else scott_bandwidth(data))
    return lambda x: kde(data, x, bandwidth=bandwidth)
