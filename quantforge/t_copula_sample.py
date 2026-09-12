"""Student-t copula sampler: dependent uniforms with heavy joint tails.

The Gaussian copula has zero tail dependence -- extremes decouple. The Student-t
copula keeps the same rank correlation but adds symmetric tail dependence that
grows as the degrees of freedom ``df`` fall, so joint crashes are far more likely.
It is the standard fix for the Gaussian copula's role in underestimating joint
default / drawdown risk.

To sample it with correlation ``R`` and ``df`` degrees of freedom:

  1. draw correlated normals ``x = L z`` (``L`` = Cholesky of ``R``);
  2. draw ``w ~ chi-square(df)`` independently;
  3. form the multivariate-t vector ``t_i = x_i * sqrt(df / w)``;
  4. map each component to a uniform with the Student-t CDF, ``u_i = T_df(t_i)``.

As ``df -> infinity`` the scaling ``sqrt(df / w) -> 1`` and this reduces to the
Gaussian copula. Uses an integer-``df`` chi-square (sum of ``df`` squared standard
normals) and a deterministic normal stream. Pure standard library.
"""

import math

from .linalg import cholesky
from .student_t import t_cdf


def _lcg_normals(seed):
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


def student_t_copula_sample(correlation, df, n, seed=1234567):
    """Draw ``n`` samples from a Student-t copula.

    Parameters
    ----------
    correlation : list[list[float]]
        Symmetric positive-definite correlation matrix (unit diagonal).
    df : int
        Degrees of freedom (>= 1). Smaller ``df`` gives heavier joint tails; large
        ``df`` approaches the Gaussian copula.
    n : int
        Number of sample vectors.
    seed : int
        Seed for the deterministic normal stream.

    Returns
    -------
    list[list[float]]
        ``n`` vectors of uniforms in (0, 1) with the target rank correlation and
        symmetric tail dependence set by ``df``.
    """
    d = len(correlation)
    if d < 1:
        raise ValueError("correlation must be non-empty")
    if any(len(row) != d for row in correlation):
        raise ValueError("correlation must be square")
    if df < 1:
        raise ValueError("df must be >= 1")
    if n < 1:
        raise ValueError("n must be >= 1")
    L = cholesky(correlation)             # raises if not positive-definite
    gen = _lcg_normals(seed)
    out = []
    for _ in range(n):
        z = [gen() for _ in range(d)]
        # Chi-square(df) as a sum of df squared standard normals.
        w = sum(gen() ** 2 for _ in range(df))
        scale = math.sqrt(df / w) if w > 0 else 1.0
        row = []
        for i in range(d):
            xi = sum(L[i][j] * z[j] for j in range(i + 1))
            row.append(t_cdf(xi * scale, df))
        out.append(row)
    return out
