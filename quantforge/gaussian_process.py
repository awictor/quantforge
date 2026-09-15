"""Gaussian process regression (exact, with the RBF and Matern kernels).

A Gaussian process places a prior over functions: any finite set of points is jointly Gaussian
with covariance given by a kernel ``k(x, x')``. Conditioning on observed ``(X, y)`` yields a
Gaussian posterior over the function value at any new input -- a mean prediction *and* a
calibrated variance, which is what makes GPs the tool of choice for Bayesian optimization and
small-data interpolation.

For ``n`` training points the exact posterior needs one Cholesky factorization of the
``n x n`` kernel matrix (plus noise on the diagonal). This module provides:

* :func:`rbf_kernel` / :func:`matern32_kernel` -- two standard stationary covariance functions.
* :func:`gp_predict` -- posterior mean and variance at test inputs.
* :func:`gp_log_marginal_likelihood` -- the evidence, for choosing kernel hyperparameters.

Inputs are scalars or equal-length coordinate lists. Pure standard library (Cholesky from
:mod:`quantforge.linalg`).
"""

import math

from .linalg import cholesky


def _sqdist(a, b):
    if isinstance(a, (list, tuple)):
        return sum((ai - bi) ** 2 for ai, bi in zip(a, b))
    return (a - b) ** 2


def rbf_kernel(a, b, length_scale=1.0, variance=1.0):
    """Squared-exponential (RBF) covariance ``variance * exp(-||a-b||^2 / (2 length_scale^2))``."""
    return variance * math.exp(-_sqdist(a, b) / (2.0 * length_scale * length_scale))


def matern32_kernel(a, b, length_scale=1.0, variance=1.0):
    """Matern-3/2 covariance ``variance (1 + sqrt3 r/l) exp(-sqrt3 r/l)`` with ``r = ||a-b||``."""
    r = math.sqrt(_sqdist(a, b))
    s = math.sqrt(3.0) * r / length_scale
    return variance * (1.0 + s) * math.exp(-s)


def _cho_solve(L, b):
    # solve (L L^T) x = b for lower-triangular L
    n = len(L)
    # forward: L z = b
    z = [0.0] * n
    for i in range(n):
        z[i] = (b[i] - sum(L[i][j] * z[j] for j in range(i))) / L[i][i]
    # backward: L^T x = z
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (z[i] - sum(L[j][i] * x[j] for j in range(i + 1, n))) / L[i][i]
    return x


def _kernel_matrix(X, kernel, length_scale, variance):
    n = len(X)
    return [[kernel(X[i], X[j], length_scale, variance) for j in range(n)] for i in range(n)]


def gp_predict(X_train, y_train, X_test, kernel=None,
               length_scale=1.0, variance=1.0, noise=1e-6):
    """Posterior mean and variance of a GP at ``X_test`` given training data ``(X_train, y_train)``.

    ``kernel`` is a covariance function ``k(a, b, length_scale, variance)``; ``noise`` is the
    observation-noise variance added to the diagonal (also regularizes the Cholesky). Returns a
    dict with ``mean`` (list) and ``var`` (list) at each test point. The prior mean is zero, so
    center ``y_train`` if it is not.
    """
    if kernel is None:
        kernel = rbf_kernel
    n = len(X_train)
    K = _kernel_matrix(X_train, kernel, length_scale, variance)
    for i in range(n):
        K[i][i] += noise
    L = cholesky(K)
    alpha = _cho_solve(L, list(y_train))
    means, variances = [], []
    for xt in X_test:
        k_star = [kernel(X_train[i], xt, length_scale, variance) for i in range(n)]
        mean = sum(k_star[i] * alpha[i] for i in range(n))
        # v = L^{-1} k_star ; var = k(xt,xt) - v.v
        v = [0.0] * n
        for i in range(n):
            v[i] = (k_star[i] - sum(L[i][j] * v[j] for j in range(i))) / L[i][i]
        kss = kernel(xt, xt, length_scale, variance)
        var = kss - sum(vi * vi for vi in v)
        means.append(mean)
        variances.append(max(var, 0.0))
    return {"mean": means, "var": variances}


def gp_log_marginal_likelihood(X_train, y_train, kernel=None,
                               length_scale=1.0, variance=1.0, noise=1e-6):
    """Log marginal likelihood ``log p(y | X)`` of the GP -- the objective for hyperparameter tuning.

    ``= -1/2 y^T K^-1 y - sum log L_ii - n/2 log(2 pi)`` with ``K = kernel + noise I = L L^T``.
    Larger is better; maximize over ``length_scale``/``variance``/``noise`` to fit the kernel.
    """
    if kernel is None:
        kernel = rbf_kernel
    n = len(X_train)
    K = _kernel_matrix(X_train, kernel, length_scale, variance)
    for i in range(n):
        K[i][i] += noise
    L = cholesky(K)
    alpha = _cho_solve(L, list(y_train))
    data_fit = -0.5 * sum(y_train[i] * alpha[i] for i in range(n))
    log_det = sum(math.log(L[i][i]) for i in range(n))   # (1/2) log det K = sum log L_ii
    return data_fit - log_det - 0.5 * n * math.log(2.0 * math.pi)
