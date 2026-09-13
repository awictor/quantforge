"""Market-stress measures: financial turbulence and the absorption ratio.

Two Kritzman-Li systemic-risk gauges:

- ``turbulence`` -- the Mahalanobis distance of a return vector from its historical
  mean given the historical covariance, ``d = (r - mu)' C^{-1} (r - mu)``. It spikes
  when returns are both large and move in unusual cross-asset combinations, so it
  flags stress that a simple volatility reading misses. Under multivariate normality
  it has mean equal to the number of assets.
- ``absorption_ratio`` -- the fraction of total variance captured by the top few
  principal components. A high ratio means risk is concentrated in a few factors
  (tightly coupled, fragile markets); a rising ratio has preceded drawdowns.

Pure standard library on top of the matrix inverse and the Jacobi eigensolver.
"""

from .portopt import _invert
from .pca import jacobi_eigen


def turbulence(observation, mean, cov):
    """Financial turbulence: Mahalanobis distance of ``observation`` from ``mean``.

    ``d = (r - mu)' C^{-1} (r - mu)`` for a single return vector ``observation`` given
    the historical ``mean`` and covariance ``cov``. Non-negative; larger means a more
    unusual (stressed) cross-asset move. Under multivariate normality ``E[d]`` equals
    the number of assets.
    """
    n = len(observation)
    if len(mean) != n or len(cov) != n or any(len(row) != n for row in cov):
        raise ValueError("dimensions of observation, mean and cov must match")
    inv = _invert(cov)
    diff = [observation[i] - mean[i] for i in range(n)]
    # d = diff' inv diff
    tmp = [sum(inv[i][j] * diff[j] for j in range(n)) for i in range(n)]
    return sum(diff[i] * tmp[i] for i in range(n))


def turbulence_series(returns, mean=None, cov=None):
    """Turbulence for each row of a return panel (rows = periods, cols = assets).

    If ``mean``/``cov`` are omitted they are estimated in-sample from ``returns``.
    Returns one turbulence value per period; the average is close to the number of
    assets when the data is multivariate normal.
    """
    m = len(returns)
    if m < 2:
        raise ValueError("need at least two observations")
    n = len(returns[0])
    if mean is None:
        mean = [sum(returns[t][i] for t in range(m)) / m for i in range(n)]
    if cov is None:
        cov = [[sum((returns[t][i] - mean[i]) * (returns[t][j] - mean[j])
                    for t in range(m)) / (m - 1) for j in range(n)] for i in range(n)]
    return [turbulence(returns[t], mean, cov) for t in range(m)]


def absorption_ratio(cov, n_factors=None):
    """Absorption ratio: variance share of the top principal components.

    Sums the ``n_factors`` largest eigenvalues of ``cov`` and divides by the total
    (the trace). ``n_factors`` defaults to about a fifth of the assets (Kritzman-Li's
    rule of thumb). In ``[0, 1]``: near ``n_factors / n`` when risk is spread evenly,
    toward 1 when a few factors dominate -- a high or rising ratio signals a fragile,
    tightly-coupled market.
    """
    n = len(cov)
    if n == 0 or any(len(row) != n for row in cov):
        raise ValueError("cov must be a non-empty square matrix")
    if n_factors is None:
        n_factors = max(1, n // 5)
    if not (1 <= n_factors <= n):
        raise ValueError("n_factors must be in [1, n]")
    eigenvalues, _ = jacobi_eigen(cov)      # sorted descending
    total = sum(max(e, 0.0) for e in eigenvalues)
    if total <= 0.0:
        raise ValueError("covariance has non-positive total variance")
    top = sum(max(eigenvalues[i], 0.0) for i in range(n_factors))
    return top / total
