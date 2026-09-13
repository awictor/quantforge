"""Compensated summation and numerically stable running statistics.

Adding many floating-point numbers loses low-order bits: a naive ``sum`` of a long or
badly-scaled series can be off by far more than machine epsilon. These recover the lost
precision:

  * ``kahan_sum`` -- Kahan compensated summation, carrying a running correction term.
  * ``neumaier_sum`` -- Neumaier's improved variant, correct even when the next term is
    larger than the running total.
  * ``accurate_dot`` -- a compensated dot product built on the same idea.
  * ``welford`` -- Welford's one-pass mean and variance, numerically stable where the
    textbook ``E[x^2] - E[x]^2`` formula catastrophically cancels for large-mean data.

Pure standard library.
"""

import math


def kahan_sum(values):
    """Kahan compensated summation of ``values``.

    Maintains a compensation term for the low-order bits dropped at each addition, so
    the result is far more accurate than the naive running sum for long or
    poorly-scaled sequences.
    """
    total = 0.0
    c = 0.0
    for v in values:
        y = v - c
        t = total + y
        c = (t - total) - y
        total = t
    return total


def neumaier_sum(values):
    """Neumaier (improved Kahan) summation.

    Like :func:`kahan_sum` but also correct when an individual term exceeds the running
    total in magnitude -- it accumulates the correction from whichever operand is
    larger. The most robust simple compensated sum.
    """
    total = 0.0
    c = 0.0
    for v in values:
        t = total + v
        if abs(total) >= abs(v):
            c += (total - t) + v
        else:
            c += (v - t) + total
        total = t
    return total + c


def accurate_dot(a, b):
    """Compensated dot product of two equal-length vectors."""
    if len(a) != len(b):
        raise ValueError("vectors must have equal length")
    total = 0.0
    c = 0.0
    for i in range(len(a)):
        p = a[i] * b[i]
        t = total + p
        if abs(total) >= abs(p):
            c += (total - t) + p
        else:
            c += (p - t) + total
        total = t
    return total + c


def welford(values):
    """Welford's stable one-pass mean and sample variance.

    Returns ``(mean, variance, n)`` with the sample (``n-1``) variance. Numerically
    stable for large-mean, small-spread data where ``mean(x^2) - mean(x)^2`` loses all
    precision to cancellation. Raises on an empty input.
    """
    n = 0
    mean = 0.0
    m2 = 0.0
    for v in values:
        n += 1
        delta = v - mean
        mean += delta / n
        m2 += delta * (v - mean)
    if n == 0:
        raise ValueError("need at least one value")
    var = m2 / (n - 1) if n > 1 else 0.0
    return mean, var, n
