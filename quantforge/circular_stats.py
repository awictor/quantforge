"""Directional (circular) statistics for angles.

Ordinary mean and variance are wrong for angles: the mean of 350 deg and 10 deg is 0 deg,
not 180 deg. Circular statistics work with the unit-vector sum instead. Provides the mean
direction, resultant length (concentration), circular variance/standard deviation, and
the Rayleigh test for non-uniformity (is there a preferred direction at all?). Angles are
in radians. Pure standard library.
"""

import math


def circular_mean(angles):
    """Mean direction of a set of angles (radians), in ``(-pi, pi]``.

    Averages the unit vectors ``(cos, sin)`` and takes the atan2 of the result, so it
    wraps correctly across the ``2*pi`` boundary. Undefined (raises) when the vectors
    cancel to the origin (no mean direction).
    """
    n = len(angles)
    if n == 0:
        raise ValueError("need at least one angle")
    c = sum(math.cos(a) for a in angles)
    s = sum(math.sin(a) for a in angles)
    if abs(c) < 1e-15 and abs(s) < 1e-15:
        raise ValueError("mean direction undefined (vectors cancel)")
    return math.atan2(s, c)


def resultant_length(angles):
    """Mean resultant length ``R`` in ``[0, 1]``: how concentrated the angles are.

    ``R = |sum e^{i theta}| / n``. ``1`` means all angles identical; ``0`` means they are
    spread so the unit vectors cancel. The basis for circular variance and the Rayleigh
    test.
    """
    n = len(angles)
    if n == 0:
        raise ValueError("need at least one angle")
    c = sum(math.cos(a) for a in angles)
    s = sum(math.sin(a) for a in angles)
    return math.hypot(c, s) / n


def circular_variance(angles):
    """Circular variance ``1 - R`` in ``[0, 1]`` (0 = concentrated, 1 = dispersed)."""
    return 1.0 - resultant_length(angles)


def circular_std(angles):
    """Circular standard deviation ``sqrt(-2 ln R)`` (radians).

    Grows without bound as the angles spread (``R -> 0``); ``0`` when all identical.
    """
    R = resultant_length(angles)
    if R <= 0.0:
        return float("inf")
    return math.sqrt(-2.0 * math.log(R))


def rayleigh_test(angles):
    """Rayleigh test for a uniform circular distribution: returns ``(R, p_value)``.

    Tests the null hypothesis that the angles are uniformly spread around the circle
    against the alternative of a single preferred direction. A small ``p_value`` rejects
    uniformity. Uses the standard ``Z = n R^2`` statistic with the Zar small-sample
    correction. Needs at least two angles.
    """
    n = len(angles)
    if n < 2:
        raise ValueError("need at least two angles")
    R = resultant_length(angles)
    Z = n * R * R
    # Zar (1999) approximation to the p-value.
    p = math.exp(-Z) * (1.0 + (2.0 * Z - Z * Z) / (4.0 * n)
                        - (24.0 * Z - 132.0 * Z ** 2 + 76.0 * Z ** 3 - 9.0 * Z ** 4)
                        / (288.0 * n * n))
    return R, min(max(p, 0.0), 1.0)
