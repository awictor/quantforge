"""Statistical power and sample-size calculations.

Given an effect size, sample size and significance level, the *power* of a test is
the probability it rejects a false null. These routines use the normal
approximation (large-sample, or the z-test), which is standard for planning:

- ``two_sample_t_power`` / ``two_sample_t_sample_size`` for a difference in means
  parameterised by Cohen's ``d = (mu1 - mu2) / sigma`` (pooled),
- ``one_sample_z_power`` / ``one_sample_z_sample_size`` for a one-sample mean shift,
- ``proportion_power`` / ``proportion_sample_size`` for a two-proportion test.

All two-sided. Pure standard library on top of the normal CDF and the inverse
error function.
"""

import math

from .mathfns import norm_cdf
from .special import erfinv


def _z(alpha_or_conf):
    """Upper-tail normal critical value for a two-sided level ``alpha``.

    ``_z(alpha)`` returns ``Phi^{-1}(1 - alpha/2)`` via ``erfinv``.
    """
    if not (0.0 < alpha_or_conf < 1.0):
        raise ValueError("level must be in (0, 1)")
    return math.sqrt(2.0) * erfinv(1.0 - alpha_or_conf)


def _norm_ppf(p):
    return math.sqrt(2.0) * erfinv(2.0 * p - 1.0)


def two_sample_t_power(effect_size, n_per_group, alpha=0.05):
    """Power of a two-sided two-sample test for a mean difference (normal approx).

    ``effect_size`` is Cohen's ``d`` (mean difference in pooled-SD units). With
    ``n_per_group`` observations in each arm the noncentrality is
    ``d sqrt(n / 2)``, and the power is ``Phi(ncp - z_{alpha/2})`` plus the far
    tail. Increases with the effect size, the sample size, and ``alpha``.
    """
    if n_per_group < 2:
        raise ValueError("n_per_group must be at least 2")
    z_a = _z(alpha)
    ncp = abs(effect_size) * math.sqrt(n_per_group / 2.0)
    return norm_cdf(ncp - z_a) + norm_cdf(-ncp - z_a)


def two_sample_t_sample_size(effect_size, power=0.80, alpha=0.05):
    """Per-group sample size for a target ``power`` in a two-sample mean test.

    Inverts the normal-approximation power: ``n = 2 (z_{alpha/2} + z_{beta})^2 /
    d^2``, rounded up. Raises for a zero effect size (infinite sample).
    """
    if effect_size == 0.0:
        raise ValueError("effect_size must be non-zero")
    if not (0.0 < power < 1.0):
        raise ValueError("power must be in (0, 1)")
    z_a = _z(alpha)
    z_b = _norm_ppf(power)
    n = 2.0 * (z_a + z_b) ** 2 / (effect_size ** 2)
    return int(math.ceil(n))


def one_sample_z_power(effect_size, n, alpha=0.05):
    """Power of a two-sided one-sample mean test (z-test, normal approx).

    ``effect_size`` is the mean shift in SD units; the noncentrality is
    ``effect_size sqrt(n)``.
    """
    if n < 1:
        raise ValueError("n must be at least 1")
    z_a = _z(alpha)
    ncp = abs(effect_size) * math.sqrt(n)
    return norm_cdf(ncp - z_a) + norm_cdf(-ncp - z_a)


def one_sample_z_sample_size(effect_size, power=0.80, alpha=0.05):
    """Sample size for a target ``power`` in a one-sample mean test.

    ``n = (z_{alpha/2} + z_{beta})^2 / d^2``, rounded up.
    """
    if effect_size == 0.0:
        raise ValueError("effect_size must be non-zero")
    if not (0.0 < power < 1.0):
        raise ValueError("power must be in (0, 1)")
    z_a = _z(alpha)
    z_b = _norm_ppf(power)
    n = (z_a + z_b) ** 2 / (effect_size ** 2)
    return int(math.ceil(n))


def proportion_power(p1, p2, n_per_group, alpha=0.05):
    """Power of a two-sided two-proportion test (normal approximation).

    Uses the unpooled standard error at the alternative and the pooled standard
    error under the null. ``p1``, ``p2`` are the two success probabilities.
    """
    if not (0.0 < p1 < 1.0 and 0.0 < p2 < 1.0):
        raise ValueError("p1 and p2 must be in (0, 1)")
    if n_per_group < 2:
        raise ValueError("n_per_group must be at least 2")
    z_a = _z(alpha)
    pbar = 0.5 * (p1 + p2)
    se0 = math.sqrt(2.0 * pbar * (1.0 - pbar) / n_per_group)
    se1 = math.sqrt(p1 * (1.0 - p1) / n_per_group + p2 * (1.0 - p2) / n_per_group)
    diff = abs(p1 - p2)
    return norm_cdf((diff - z_a * se0) / se1)


def proportion_sample_size(p1, p2, power=0.80, alpha=0.05):
    """Per-group sample size for a target ``power`` in a two-proportion test.

    ``n = (z_{alpha/2} sqrt(2 pbar (1-pbar)) + z_{beta} sqrt(p1(1-p1)+p2(1-p2)))^2 /
    (p1 - p2)^2``, rounded up.
    """
    if not (0.0 < p1 < 1.0 and 0.0 < p2 < 1.0):
        raise ValueError("p1 and p2 must be in (0, 1)")
    if p1 == p2:
        raise ValueError("p1 and p2 must differ")
    if not (0.0 < power < 1.0):
        raise ValueError("power must be in (0, 1)")
    z_a = _z(alpha)
    z_b = _norm_ppf(power)
    pbar = 0.5 * (p1 + p2)
    term = (z_a * math.sqrt(2.0 * pbar * (1.0 - pbar))
            + z_b * math.sqrt(p1 * (1.0 - p1) + p2 * (1.0 - p2)))
    n = term ** 2 / (p1 - p2) ** 2
    return int(math.ceil(n))
