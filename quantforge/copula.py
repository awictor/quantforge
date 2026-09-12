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


def vasicek_loss_cdf(loss, pd, rho):
    """CDF of the large-homogeneous-portfolio loss fraction (Vasicek limit).

    In the single-factor Gaussian-copula limit of an infinitely granular pool with
    default probability ``pd`` and asset correlation ``rho``, the fractional loss
    ``L`` has closed-form CDF

        P(L <= x) = Phi( (sqrt(1 - rho) Phi^{-1}(x) - Phi^{-1}(pd)) / sqrt(rho) ).

    ``loss`` is a fraction in ``[0, 1]`` (LGD assumed 1). Increasing in ``loss``.
    """
    if not (0.0 <= loss <= 1.0):
        raise ValueError("loss must be in [0, 1]")
    if not (0.0 < pd < 1.0):
        raise ValueError("pd must be in (0, 1)")
    if not (0.0 < rho < 1.0):
        raise ValueError("rho must be in (0, 1)")
    if loss <= 0.0:
        return 0.0
    if loss >= 1.0:
        return 1.0
    num = math.sqrt(1.0 - rho) * norm_ppf(loss) - norm_ppf(pd)
    return norm_cdf(num / math.sqrt(rho))


def vasicek_loss_quantile(q, pd, rho):
    """Portfolio loss at confidence ``q`` (the Vasicek/Basel capital formula).

    Inverse of :func:`vasicek_loss_cdf`:

        L(q) = Phi( (Phi^{-1}(pd) + sqrt(rho) Phi^{-1}(q)) / sqrt(1 - rho) ).

    The worst-case loss not exceeded with probability ``q`` -- the basis of the
    Basel IRB capital charge. Increasing in ``q``, ``pd`` and ``rho``.
    """
    if not (0.0 < q < 1.0):
        raise ValueError("q must be in (0, 1)")
    if not (0.0 < pd < 1.0):
        raise ValueError("pd must be in (0, 1)")
    if not (0.0 < rho < 1.0):
        raise ValueError("rho must be in (0, 1)")
    num = norm_ppf(pd) + math.sqrt(rho) * norm_ppf(q)
    return norm_cdf(num / math.sqrt(1.0 - rho))


def cdo_tranche_expected_loss(attachment, detachment, pd, rho, n_steps=2000):
    """Expected loss of a CDO tranche in the Vasicek large-pool limit.

    Integrates the portfolio loss distribution over the tranche
    ``[attachment, detachment]`` and normalizes by the tranche width, giving the
    expected tranche loss as a fraction of the tranche notional. Equity (low
    attachment) tranches lose more than senior tranches at the same correlation.
    Trapezoidal integration of ``E[min(max(L - a, 0), d - a)] / (d - a)`` using
    the survival ``1 - F(l)``.
    """
    if not (0.0 <= attachment < detachment <= 1.0):
        raise ValueError("require 0 <= attachment < detachment <= 1")
    width = detachment - attachment
    # E[tranche loss] = integral_a^d (1 - F(l)) dl  (layer expected loss).
    dl = (detachment - attachment) / n_steps
    total = 0.0
    for k in range(n_steps + 1):
        l = attachment + k * dl
        surv = 1.0 - vasicek_loss_cdf(l, pd, rho)
        w = 0.5 if (k == 0 or k == n_steps) else 1.0
        total += w * surv * dl
    return total / width


def cdo_tranche_expected_loss_mc(attachment, detachment, pd, rho, n_names=100,
                                 n_paths=20000, seed=8675309):
    """Monte Carlo CDO tranche expected loss under the single-factor model.

    Simulates a finite pool of ``n_names``: a common factor ``M`` and idiosyncratic
    shocks give each name's asset value ``sqrt(rho) M + sqrt(1 - rho) Z_i``; a name
    defaults when it falls below ``Phi^{-1}(pd)``. Averages the tranche loss over
    the portfolio-loss fraction across paths. An independent finite-pool reference
    for the large-pool :func:`cdo_tranche_expected_loss` (they agree as
    ``n_names -> inf``). Deterministic per seed.
    """
    if not (0.0 <= attachment < detachment <= 1.0):
        raise ValueError("require 0 <= attachment < detachment <= 1")
    if not (0.0 < pd < 1.0):
        raise ValueError("pd must be in (0, 1)")
    if not (0.0 < rho < 1.0):
        raise ValueError("rho must be in (0, 1)")
    if n_names < 1 or n_paths < 1:
        raise ValueError("n_names and n_paths must be positive")
    threshold = norm_ppf(pd)
    sq_rho = math.sqrt(rho)
    sq_1mrho = math.sqrt(1.0 - rho)
    width = detachment - attachment
    state = seed & 0xFFFFFFFF

    def _unif():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return (state + 0.5) / 0x80000000

    def _normal():
        u1 = _unif()
        u2 = _unif()
        return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)

    total = 0.0
    for _ in range(n_paths):
        m = _normal()
        defaults = 0
        for _i in range(n_names):
            asset = sq_rho * m + sq_1mrho * _normal()
            if asset < threshold:
                defaults += 1
        loss = defaults / n_names
        tranche_loss = min(max(loss - attachment, 0.0), width)
        total += tranche_loss
    return total / n_paths / width


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
