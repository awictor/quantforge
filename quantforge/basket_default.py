"""kth-to-default basket probabilities (one-factor Gaussian copula).

In a homogeneous basket of ``n`` names, each with default probability ``pd`` and
pairwise (asset-correlation) ``rho``, defaults are driven by a common factor ``M``:
name ``i`` defaults when its latent variable ``sqrt(rho) M + sqrt(1-rho) Z_i``
falls below the default threshold. Conditional on ``M`` the names are independent,
so the number of defaults is binomial with a factor-dependent probability; the
unconditional distribution is that binomial integrated over the standard-normal
``M`` (Gauss-Hermite / midpoint quadrature).

``basket_default_distribution`` returns ``P(exactly k defaults)`` for
``k = 0..n``; ``kth_to_default_probability`` gives ``P(at least k defaults)`` --
the payoff trigger of a kth-to-default swap. Pure standard library.
"""

import math

from .mathfns import norm_cdf, norm_ppf


def _cond_pd(pd, rho, m):
    """Default probability of a name conditional on the common factor ``m``."""
    if pd <= 0.0:
        return 0.0
    if pd >= 1.0:
        return 1.0
    c = norm_ppf(pd)
    return norm_cdf((c - math.sqrt(rho) * m) / math.sqrt(1.0 - rho))


def basket_default_distribution(n, pd, rho, n_quad=200):
    """Distribution of the number of defaults in a homogeneous basket.

    Returns a list ``p`` of length ``n + 1`` with ``p[k] = P(exactly k defaults)``,
    integrating the conditional binomial over the common factor by midpoint
    quadrature on ``[-8, 8]``. The list sums to 1.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if not (0.0 <= pd <= 1.0):
        raise ValueError("pd must be in [0, 1]")
    if not (0.0 <= rho < 1.0):
        raise ValueError("rho must be in [0, 1)")

    lo, hi = -8.0, 8.0
    dm = (hi - lo) / n_quad
    dist = [0.0] * (n + 1)
    for i in range(n_quad):
        m = lo + (i + 0.5) * dm
        phi = math.exp(-0.5 * m * m) / math.sqrt(2.0 * math.pi)
        w = phi * dm
        q = _cond_pd(pd, rho, m)
        # Conditional binomial pmf via the recurrence.
        pmf = (1.0 - q) ** n if q < 1.0 else 0.0
        ratio = (q / (1.0 - q)) if 0.0 < q < 1.0 else None
        if q >= 1.0:
            dist[n] += w
            continue
        dist[0] += w * pmf
        for k in range(1, n + 1):
            pmf *= ratio * (n - k + 1) / k
            dist[k] += w * pmf
    total = sum(dist)
    return [d / total for d in dist]


def kth_to_default_probability(n, k, pd, rho, n_quad=200):
    """Probability of at least ``k`` defaults in a homogeneous basket.

    The trigger probability of a kth-to-default swap. Monotone: decreasing in
    ``k`` (harder to reach more defaults) and, for a fixed ``k > 1``, increasing
    in correlation (defaults cluster). ``k = 1`` is first-to-default.
    """
    if not (1 <= k <= n):
        raise ValueError("require 1 <= k <= n")
    dist = basket_default_distribution(n, pd, rho, n_quad)
    return sum(dist[k:])
