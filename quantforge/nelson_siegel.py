"""Nelson-Siegel and Svensson parametric yield curves.

The Nelson-Siegel (1987) zero-rate curve fits a level, slope, and curvature with
one decay parameter; the Svensson (1994) extension adds a second curvature hump.
Both are the workhorse parametric forms central banks publish. This module
evaluates the zero rate, discount factor, and instantaneous forward rate. Pure
standard library.
"""

import math


def nelson_siegel_zero(t, beta0, beta1, beta2, tau):
    """Nelson-Siegel zero rate at maturity ``t``.

    ``z(t) = beta0 + (beta1 + beta2) (1 - e^{-t/tau}) / (t/tau) - beta2 e^{-t/tau}``.
    ``beta0`` is the long-run level, ``beta0 + beta1`` the short rate (``t -> 0``),
    and ``beta2`` scales the medium-term curvature hump with decay ``tau``.
    """
    if tau <= 0:
        raise ValueError("tau must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if t == 0.0:
        return beta0 + beta1
    x = t / tau
    loading = (1.0 - math.exp(-x)) / x
    return beta0 + beta1 * loading + beta2 * (loading - math.exp(-x))


def svensson_zero(t, beta0, beta1, beta2, beta3, tau1, tau2):
    """Svensson zero rate: Nelson-Siegel plus a second curvature term.

    Adds ``beta3 ((1 - e^{-t/tau2})/(t/tau2) - e^{-t/tau2})`` with its own decay
    ``tau2`` for a second hump. Reduces to :func:`nelson_siegel_zero` when
    ``beta3 = 0``.
    """
    if tau1 <= 0 or tau2 <= 0:
        raise ValueError("tau1 and tau2 must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    base = nelson_siegel_zero(t, beta0, beta1, beta2, tau1)
    if t == 0.0:
        return base
    x2 = t / tau2
    loading2 = (1.0 - math.exp(-x2)) / x2
    return base + beta3 * (loading2 - math.exp(-x2))


def nelson_siegel_discount(t, beta0, beta1, beta2, tau):
    """Discount factor ``exp(-z(t) t)`` from the Nelson-Siegel zero rate."""
    if t <= 0:
        return 1.0
    return math.exp(-nelson_siegel_zero(t, beta0, beta1, beta2, tau) * t)


def nelson_siegel_forward(t, beta0, beta1, beta2, tau):
    """Instantaneous forward rate under Nelson-Siegel.

    ``f(t) = beta0 + beta1 e^{-t/tau} + beta2 (t/tau) e^{-t/tau}``. Equals the
    short rate ``beta0 + beta1`` at ``t = 0`` and the long level ``beta0`` as
    ``t -> inf``.
    """
    if tau <= 0:
        raise ValueError("tau must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    x = t / tau
    e = math.exp(-x)
    return beta0 + beta1 * e + beta2 * x * e
