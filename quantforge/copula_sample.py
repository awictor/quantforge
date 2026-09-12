"""Gaussian-copula sampler: draw correlated uniforms with a target correlation.

A Gaussian copula couples arbitrary marginals through a multivariate-normal
dependence structure. To sample it:

  1. draw independent standard normals ``z``;
  2. correlate them via the Cholesky factor ``L`` of the target correlation matrix
     ``R`` -> ``x = L z`` has correlation ``R``;
  3. map each component back to a uniform with the normal CDF, ``u = Phi(x)``.

The resulting ``u`` vectors are uniform on each margin but dependent, with a rank
correlation close to ``R``. Feeding each margin through an inverse-CDF then yields
correlated draws from any target distributions -- the standard way to simulate a
dependent portfolio. A deterministic linear-congruential stream keeps results
reproducible per seed. Pure standard library.
"""

import math

from .linalg import cholesky
from .mathfns import norm_ppf, norm_cdf


def _lcg_normals(seed):
    """Deterministic standard-normal generator (LCG uniforms -> Box-Muller)."""
    state = seed & 0x7FFFFFFF
    cache = []

    def _uni():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return (state + 0.5) / 0x80000000

    def _next():
        if cache:
            return cache.pop()
        u1 = _uni()
        u2 = _uni()
        r = math.sqrt(-2.0 * math.log(u1))
        cache.append(r * math.sin(2.0 * math.pi * u2))
        return r * math.cos(2.0 * math.pi * u2)

    return _next


def gaussian_copula_sample(correlation, n, seed=1234567):
    """Draw ``n`` samples from a Gaussian copula with the given correlation matrix.

    Parameters
    ----------
    correlation : list[list[float]]
        A symmetric positive-definite correlation matrix (unit diagonal).
    n : int
        Number of sample vectors to draw.
    seed : int
        Seed for the deterministic normal stream.

    Returns
    -------
    list[list[float]]
        ``n`` vectors of uniforms in (0, 1); each margin is uniform and the
        cross-margin rank correlation approximates ``correlation``.
    """
    d = len(correlation)
    if d < 1:
        raise ValueError("correlation must be non-empty")
    if any(len(row) != d for row in correlation):
        raise ValueError("correlation must be square")
    if n < 1:
        raise ValueError("n must be >= 1")
    L = cholesky(correlation)             # raises if not positive-definite
    gen = _lcg_normals(seed)
    out = []
    for _ in range(n):
        z = [gen() for _ in range(d)]
        # x = L z, then map to uniform via the normal CDF.
        row = []
        for i in range(d):
            xi = sum(L[i][j] * z[j] for j in range(i + 1))
            row.append(norm_cdf(xi))
        out.append(row)
    return out


def inverse_transform(u, ppf):
    """Map copula uniforms to a target margin via its inverse CDF ``ppf``.

    ``u`` is a sequence of uniforms in (0, 1); ``ppf`` maps a probability to a
    quantile (e.g. :func:`quantforge.norm_ppf` for a normal margin). Returns the
    transformed sample.
    """
    return [ppf(ui) for ui in u]
