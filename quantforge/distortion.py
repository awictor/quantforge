"""Distortion (spectral) risk pricing: Wang transform and proportional hazard.

A distortion risk measure prices a loss ``S`` by integrating its distorted
survival function,

    premium = integral_0^inf g(P(S > x)) dx,

where ``g: [0,1] -> [0,1]`` is increasing with ``g(0)=0, g(1)=1``. A concave ``g``
loads the tail, giving a premium above the expected loss. Two standard distortions:

  * Wang transform: ``g(u) = Phi(Phi^{-1}(u) + lambda)`` -- shifts the loss
    distribution by ``lambda`` standard-normal units; ``lambda = 0`` recovers the
    mean, positive ``lambda`` adds a risk load. Preserves normality.
  * Proportional hazard: ``g(u) = u^{1/rho}`` for ``rho >= 1`` -- inflates the
    survival probabilities, loading the tail more as ``rho`` grows.

Works on a discrete loss distribution (probabilities on an integer grid). Pure
standard library.
"""

from .mathfns import norm_cdf, norm_ppf


def _survival(pmf):
    """Survival ``S_k = P(loss > k)`` on the grid from a pmf list."""
    n = len(pmf)
    surv = [0.0] * n
    running = 0.0
    for k in range(n - 1, -1, -1):
        surv[k] = running          # P(loss > k) = sum_{j>k} pmf[j]
        running += pmf[k]
    return surv


def wang_premium(pmf, lam):
    """Wang-transform premium of a discrete loss distribution.

    Distorts the survival function by ``g(u) = Phi(Phi^{-1}(u) + lam)`` and sums
    it over the grid (unit spacing). ``lam = 0`` gives the expected loss; positive
    ``lam`` adds a risk load. Monotone increasing in ``lam``.
    """
    surv = _survival(pmf)
    total = 0.0
    for s in surv:
        if s <= 0.0:
            continue
        if s >= 1.0:
            total += 1.0
            continue
        total += norm_cdf(norm_ppf(s) + lam)
    return total


def proportional_hazard_premium(pmf, rho):
    """Proportional-hazard distortion premium ``sum_k S_k^{1/rho}``.

    ``rho = 1`` gives the expected loss; ``rho > 1`` loads the tail. Requires
    ``rho >= 1``.
    """
    if rho < 1.0:
        raise ValueError("rho must be >= 1")
    surv = _survival(pmf)
    return sum(s ** (1.0 / rho) for s in surv if s > 0.0)


def expected_loss(pmf):
    """Expected loss ``sum_k k * pmf[k]`` (= sum of the survival function)."""
    return sum(k * pmf[k] for k in range(len(pmf)))
