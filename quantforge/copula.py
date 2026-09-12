"""Bivariate copulas: Gaussian and Archimedean (Clayton, Gumbel).

A copula ``C(u, v)`` is a joint CDF on the unit square with uniform margins,
capturing the dependence structure separately from the marginals. This module
provides the Gaussian copula and two Archimedean families, their tail-dependence
coefficients, and the Kendall's-tau relationships used to calibrate them. Every
copula satisfies the boundary conditions ``C(u, 0) = 0``, ``C(u, 1) = u`` and
reduces to the independence copula ``u v`` at zero dependence. Pure standard
library.
"""

import math

from .mathfns import norm_cdf, norm_ppf


def _bivariate_normal_cdf(x, y, rho):
    """Standard bivariate normal CDF ``P(X <= x, Y <= y)`` via Drezner-Wesolowsky.

    Gauss-Legendre quadrature of the standard formula; accurate to ~1e-10 for the
    dependence work here. Handles the ``|rho| -> 1`` limits by comonotone/
    countermonotone bounds.
    """
    if rho <= -1.0:
        return max(norm_cdf(x) + norm_cdf(y) - 1.0, 0.0)
    if rho >= 1.0:
        return min(norm_cdf(x), norm_cdf(y))
    if x == 0.0 and y == 0.0:
        return 0.25 + math.asin(rho) / (2.0 * math.pi)
    # Drezner-Wesolowsky 5-point Gauss quadrature on the reduced integral.
    nodes = [0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992]
    weights = [0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042]
    hk = x * y
    bvn = 0.0
    for w, t in zip(weights, nodes):
        r = rho * t
        denom = 1.0 - r * r
        bvn += w * math.exp((r * hk - 0.5 * (x * x + y * y)) / denom) / math.sqrt(denom)
    bvn = bvn * rho + norm_cdf(x) * norm_cdf(y)
    return min(max(bvn, 0.0), 1.0)


def gaussian_copula(u, v, rho):
    """Gaussian copula ``C(u, v) = Phi_rho(Phi^{-1}(u), Phi^{-1}(v))``.

    The dependence structure of a bivariate normal with correlation ``rho``. Zero
    tail dependence for ``|rho| < 1``; reduces to ``u v`` at ``rho = 0``.
    """
    if not (0.0 <= u <= 1.0 and 0.0 <= v <= 1.0):
        raise ValueError("u and v must be in [0, 1]")
    if u == 0.0 or v == 0.0:
        return 0.0
    if u == 1.0:
        return v
    if v == 1.0:
        return u
    return _bivariate_normal_cdf(norm_ppf(u), norm_ppf(v), rho)


def clayton_copula(u, v, theta):
    """Clayton copula ``(u^{-theta} + v^{-theta} - 1)^{-1/theta}`` (``theta > 0``).

    Lower-tail dependent (assets crash together); reduces to independence as
    ``theta -> 0``.
    """
    if not (0.0 <= u <= 1.0 and 0.0 <= v <= 1.0):
        raise ValueError("u and v must be in [0, 1]")
    if theta <= 0.0:
        raise ValueError("theta must be positive")
    if u == 0.0 or v == 0.0:
        return 0.0
    if u == 1.0:
        return v
    if v == 1.0:
        return u
    return (u ** (-theta) + v ** (-theta) - 1.0) ** (-1.0 / theta)


def gumbel_copula(u, v, theta):
    """Gumbel copula ``exp(-((-ln u)^theta + (-ln v)^theta)^{1/theta})`` (``theta >= 1``).

    Upper-tail dependent (assets rally together); reduces to independence at
    ``theta = 1``.
    """
    if not (0.0 <= u <= 1.0 and 0.0 <= v <= 1.0):
        raise ValueError("u and v must be in [0, 1]")
    if theta < 1.0:
        raise ValueError("theta must be >= 1")
    if u == 0.0 or v == 0.0:
        return 0.0
    if u == 1.0:
        return v
    if v == 1.0:
        return u
    lu = (-math.log(u)) ** theta
    lv = (-math.log(v)) ** theta
    return math.exp(-(lu + lv) ** (1.0 / theta))


def frank_copula(u, v, theta):
    """Frank copula (``theta != 0``), symmetric with no tail dependence.

    ``C(u, v) = -1/theta * ln(1 + (e^{-theta u} - 1)(e^{-theta v} - 1) /
    (e^{-theta} - 1))``. Positive dependence for ``theta > 0``, negative for
    ``theta < 0``; reduces to independence as ``theta -> 0``.
    """
    if not (0.0 <= u <= 1.0 and 0.0 <= v <= 1.0):
        raise ValueError("u and v must be in [0, 1]")
    if u == 0.0 or v == 0.0:
        return 0.0
    if u == 1.0:
        return v
    if v == 1.0:
        return u
    if abs(theta) < 1e-8:
        return u * v
    num = (math.exp(-theta * u) - 1.0) * (math.exp(-theta * v) - 1.0)
    den = math.exp(-theta) - 1.0
    return -1.0 / theta * math.log1p(num / den)


def gaussian_copula_joint_default(pd1, pd2, rho):
    """Joint default probability of two names under the Gaussian copula.

    Both default when their latent normals fall below their default thresholds
    ``Phi^{-1}(pd_i)``; the joint probability is the Gaussian copula
    ``C(pd1, pd2; rho)``. Rises above the independent product ``pd1 * pd2`` for
    ``rho > 0`` and equals it at ``rho = 0``.
    """
    if not (0.0 <= pd1 <= 1.0 and 0.0 <= pd2 <= 1.0):
        raise ValueError("default probabilities must be in [0, 1]")
    return gaussian_copula(pd1, pd2, rho)


def first_to_default_probability(pd1, pd2, rho):
    """Probability that at least one of two names defaults (Gaussian copula).

    ``P(A or B) = pd1 + pd2 - C(pd1, pd2; rho)`` by inclusion-exclusion. Lies
    between ``max(pd1, pd2)`` and ``min(pd1 + pd2, 1)``, and falls as correlation
    rises (correlated defaults overlap more, so fewer *distinct* default events).
    """
    joint = gaussian_copula_joint_default(pd1, pd2, rho)
    return pd1 + pd2 - joint


def clayton_lower_tail_dependence(theta):
    """Lower-tail dependence of the Clayton copula ``2^{-1/theta}``.

    In ``(0, 1)`` for ``theta > 0`` -- rising toward 1 as ``theta`` grows (stronger
    joint-crash dependence).
    """
    if theta <= 0.0:
        raise ValueError("theta must be positive")
    return 2.0 ** (-1.0 / theta)


def gumbel_upper_tail_dependence(theta):
    """Upper-tail dependence of the Gumbel copula ``2 - 2^{1/theta}``.

    Zero at ``theta = 1`` (independence) rising toward 1 as ``theta -> inf``.
    """
    if theta < 1.0:
        raise ValueError("theta must be >= 1")
    return 2.0 - 2.0 ** (1.0 / theta)


def clayton_theta_from_tau(tau):
    """Clayton ``theta`` from Kendall's tau: ``2 tau / (1 - tau)``.

    Inverse of ``tau = theta / (theta + 2)``. Requires ``0 <= tau < 1``.
    """
    if not (0.0 <= tau < 1.0):
        raise ValueError("tau must be in [0, 1)")
    return 2.0 * tau / (1.0 - tau)


def gumbel_theta_from_tau(tau):
    """Gumbel ``theta`` from Kendall's tau: ``1 / (1 - tau)``.

    Inverse of ``tau = 1 - 1/theta``. Requires ``0 <= tau < 1``.
    """
    if not (0.0 <= tau < 1.0):
        raise ValueError("tau must be in [0, 1)")
    return 1.0 / (1.0 - tau)
