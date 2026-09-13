"""Multivariate normal density and log-determinant via Cholesky.

The multivariate-normal log-density needs the log-determinant and the quadratic form
of the covariance's inverse. Both are computed stably from the Cholesky factor
``Sigma = L L'``:

- ``log_determinant(Sigma) = 2 sum log L_ii`` (no overflow, unlike a raw product),
- the Mahalanobis term ``(x - mu)' Sigma^{-1} (x - mu) = ||L^{-1}(x - mu)||^2`` by a
  triangular solve.

``mvn_logpdf`` assembles the full log-density. Pure standard library on top of the
Cholesky factorization.
"""

import math

from .linalg import cholesky


def log_determinant(cov):
    """Log-determinant of a symmetric positive-definite matrix via Cholesky.

    ``log det Sigma = 2 sum_i log L_ii``. Numerically stable where a direct product
    of eigenvalues (or the determinant) would under/overflow. Raises if ``cov`` is
    not positive definite.
    """
    L = cholesky(cov)
    return 2.0 * sum(math.log(L[i][i]) for i in range(len(L)))


def _forward_solve(L, b):
    """Solve L y = b for lower-triangular L."""
    n = len(L)
    y = [0.0] * n
    for i in range(n):
        y[i] = (b[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]
    return y


def mvn_logpdf(x, mean, cov):
    """Log-density of the multivariate normal ``N(mean, cov)`` at ``x``.

    ``-0.5 [ k ln(2 pi) + ln|Sigma| + (x-mu)' Sigma^{-1} (x-mu) ]``. The quadratic
    form is evaluated as ``||L^{-1}(x-mu)||^2`` from the Cholesky factor, avoiding an
    explicit inverse. Reduces to the univariate normal log-density for ``k = 1``.
    """
    k = len(x)
    if len(mean) != k or len(cov) != k or any(len(row) != k for row in cov):
        raise ValueError("dimensions of x, mean and cov must match")
    L = cholesky(cov)
    diff = [x[i] - mean[i] for i in range(k)]
    y = _forward_solve(L, diff)                    # L y = diff
    quad = sum(v * v for v in y)                   # (x-mu)' Sigma^-1 (x-mu)
    logdet = 2.0 * sum(math.log(L[i][i]) for i in range(k))
    return -0.5 * (k * math.log(2.0 * math.pi) + logdet + quad)


def mvn_pdf(x, mean, cov):
    """Density of the multivariate normal ``N(mean, cov)`` at ``x`` (``exp`` of the log)."""
    return math.exp(mvn_logpdf(x, mean, cov))
