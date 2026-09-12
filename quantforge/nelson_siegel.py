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


def _ns_loadings(t, tau):
    """The three Nelson-Siegel factor loadings at maturity ``t`` for a fixed tau."""
    if t == 0.0:
        return 1.0, 1.0, 0.0
    x = t / tau
    load = (1.0 - math.exp(-x)) / x
    return 1.0, load, load - math.exp(-x)


def fit_nelson_siegel(maturities, zero_rates, tau_grid=None):
    """Least-squares fit of Nelson-Siegel parameters to observed zero rates.

    For a fixed decay ``tau`` the three betas enter linearly (the level/slope/
    curvature loadings), so they are solved by ordinary least squares; ``tau`` is
    chosen by a grid search minimizing the residual sum of squares. Returns
    ``(beta0, beta1, beta2, tau)``. Recovers the true parameters exactly on
    noiseless data whose ``tau`` is in the grid.
    """
    from .factor_model import factor_regression
    if len(maturities) != len(zero_rates):
        raise ValueError("maturities and zero_rates must have equal length")
    if len(maturities) < 3:
        raise ValueError("need at least three points to fit three betas")
    if tau_grid is None:
        tau_grid = [0.25 * k for k in range(1, 41)]   # 0.25 .. 10 years
    best = None
    for tau in tau_grid:
        # Loadings 2 and 3 as regressors; loading 1 is the intercept (all ones).
        f2 = [_ns_loadings(t, tau)[1] for t in maturities]
        f3 = [_ns_loadings(t, tau)[2] for t in maturities]
        fit = factor_regression(list(zero_rates), [f2, f3])
        beta0 = fit["alpha"]
        beta1, beta2 = fit["betas"]
        sse = 0.0
        for t, z in zip(maturities, zero_rates):
            sse += (nelson_siegel_zero(t, beta0, beta1, beta2, tau) - z) ** 2
        if best is None or sse < best[0]:
            best = (sse, beta0, beta1, beta2, tau)
    return best[1], best[2], best[3], best[4]


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
