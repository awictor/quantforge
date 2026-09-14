"""Von Mises distribution: the circular analogue of the normal.

The von Mises distribution is the natural bell curve on a circle -- the maximum-entropy
distribution for a given mean direction and concentration. Its density is
``exp(kappa cos(theta - mu)) / (2 pi I0(kappa))`` where ``mu`` is the mean angle and
``kappa`` the concentration (``kappa = 0`` is uniform, large ``kappa`` is a tight peak).
Provides the density, the modified Bessel functions it needs, and a maximum-likelihood
fit of ``(mu, kappa)`` from data. Pure standard library.
"""

import math

from .circular_stats import circular_mean, resultant_length


def bessel_i0(x):
    """Modified Bessel function of the first kind, order 0, ``I0(x)``.

    Uses the Abramowitz & Stegun polynomial approximations (accurate to ~1e-7), the same
    ones used for the von Mises normalizing constant.
    """
    ax = abs(x)
    if ax < 3.75:
        t = (x / 3.75) ** 2
        return (1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (
            0.2659732 + t * (0.0360768 + t * 0.0045813))))))
    t = 3.75 / ax
    return (math.exp(ax) / math.sqrt(ax)) * (0.39894228 + t * (0.01328592 + t * (
        0.00225319 + t * (-0.00157565 + t * (0.00916281 + t * (-0.02057706 + t * (
            0.02635537 + t * (-0.01647633 + t * 0.00392377))))))))


def bessel_i1(x):
    """Modified Bessel function of the first kind, order 1, ``I1(x)`` (A&S approximation)."""
    ax = abs(x)
    if ax < 3.75:
        t = (x / 3.75) ** 2
        val = ax * (0.5 + t * (0.87890594 + t * (0.51498869 + t * (0.15084934 + t * (
            0.02658733 + t * (0.00301532 + t * 0.00032411))))))
    else:
        t = 3.75 / ax
        val = (math.exp(ax) / math.sqrt(ax)) * (0.39894228 + t * (-0.03988024 + t * (
            -0.00362018 + t * (0.00163801 + t * (-0.01031555 + t * (0.02282967 + t * (
                -0.02895312 + t * (0.01787654 + t * -0.00420059))))))))
    return -val if x < 0 else val


def von_mises_pdf(theta, mu, kappa):
    """Von Mises density at angle ``theta`` with mean ``mu`` and concentration ``kappa``.

    ``exp(kappa cos(theta - mu)) / (2 pi I0(kappa))``. ``kappa >= 0``; ``kappa = 0`` gives
    the uniform density ``1 / (2 pi)``. Integrates to 1 over any ``2*pi`` interval.
    """
    if kappa < 0:
        raise ValueError("kappa must be non-negative")
    return math.exp(kappa * math.cos(theta - mu)) / (2.0 * math.pi * bessel_i0(kappa))


def von_mises_fit(angles, tol=1e-10, max_iter=100):
    """Maximum-likelihood fit of a von Mises distribution to ``angles`` (radians).

    The MLE mean direction is the sample :func:`circular_mean`; the concentration
    ``kappa`` solves ``I1(kappa)/I0(kappa) = R`` (the mean resultant length), found here by
    Newton's method with a standard closed-form seed. Returns ``(mu, kappa)``.
    """
    if len(angles) < 2:
        raise ValueError("need at least two angles")
    mu = circular_mean(angles)
    R = resultant_length(angles)
    if R <= 0.0:
        return mu, 0.0
    if R >= 1.0:
        return mu, float("inf")
    # Closed-form initial estimate (Fisher 1993).
    if R < 0.53:
        kappa = 2.0 * R + R ** 3 + 5.0 * R ** 5 / 6.0
    elif R < 0.85:
        kappa = -0.4 + 1.39 * R + 0.43 / (1.0 - R)
    else:
        kappa = 1.0 / (R ** 3 - 4.0 * R ** 2 + 3.0 * R)
    # Newton refine on A(kappa) = I1/I0 - R.
    for _ in range(max_iter):
        i0 = bessel_i0(kappa)
        i1 = bessel_i1(kappa)
        A = i1 / i0
        # A'(kappa) = 1 - A^2 - A/kappa
        Aprime = 1.0 - A * A - A / kappa
        if abs(Aprime) < 1e-15:
            break
        step = (A - R) / Aprime
        kappa -= step
        if kappa <= 0:
            kappa = 1e-6
        if abs(step) < tol:
            break
    return mu, kappa
