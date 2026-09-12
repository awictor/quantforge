"""Ledoit-Wolf shrinkage covariance estimation.

The sample covariance matrix is noisy when the number of assets ``p`` is not small
relative to the number of observations ``n`` -- its extreme eigenvalues are biased
and it may be ill-conditioned or singular. Ledoit and Wolf (2004) shrink it toward
a structured target ``F`` (here the constant-correlation matrix) by a data-driven
intensity ``delta`` in ``[0, 1]``:

    Sigma_hat = delta * F + (1 - delta) * S

where ``S`` is the sample covariance. The optimal ``delta`` minimizes the expected
Frobenius distance to the true covariance and has the closed form
``delta* = (pi - rho) / gamma / n`` (clipped to ``[0, 1]``), with ``pi`` the sum of
asymptotic variances of the sample-covariance entries, ``rho`` their covariance
with the target, and ``gamma`` the squared Frobenius misfit of ``F`` to ``S``.

The result is always at least as well-conditioned as ``S`` and is positive
semidefinite when the target is. Pure standard library.
"""

import math


def _demean(returns):
    """Return (columns, n, p, means): column-major data and per-column means."""
    n = len(returns)
    if n < 2:
        raise ValueError("need at least 2 observations")
    p = len(returns[0])
    if any(len(row) != p for row in returns):
        raise ValueError("all observations must have the same length")
    means = [sum(returns[t][i] for t in range(n)) / n for i in range(p)]
    return n, p, means


def sample_covariance(returns):
    """Maximum-likelihood sample covariance (divisor ``n``) of a return matrix.

    ``returns`` is a sequence of ``n`` observations, each a length-``p`` sequence.
    Returns a ``p x p`` list-of-lists.
    """
    n, p, means = _demean(returns)
    cov = [[0.0] * p for _ in range(p)]
    for t in range(n):
        dev = [returns[t][i] - means[i] for i in range(p)]
        for i in range(p):
            di = dev[i]
            row = cov[i]
            for j in range(i, p):
                row[j] += di * dev[j]
    for i in range(p):
        for j in range(i, p):
            cov[i][j] /= n
            cov[j][i] = cov[i][j]
    return cov


def constant_correlation_target(cov):
    """Constant-correlation shrinkage target from a covariance matrix.

    Keeps each asset's own variance but replaces every pairwise correlation with
    the average sample correlation ``r_bar``: ``F_ij = r_bar * sqrt(S_ii S_jj)``
    for ``i != j`` and ``F_ii = S_ii``.
    """
    p = len(cov)
    std = [math.sqrt(cov[i][i]) for i in range(p)]
    # Average off-diagonal correlation.
    s = 0.0
    cnt = 0
    for i in range(p):
        for j in range(i + 1, p):
            denom = std[i] * std[j]
            if denom > 0.0:
                s += cov[i][j] / denom
                cnt += 1
    r_bar = s / cnt if cnt else 0.0
    F = [[0.0] * p for _ in range(p)]
    for i in range(p):
        F[i][i] = cov[i][i]
        for j in range(i + 1, p):
            F[i][j] = F[j][i] = r_bar * std[i] * std[j]
    return F, r_bar


def ledoit_wolf_shrinkage(returns):
    """Ledoit-Wolf (2004) constant-correlation shrinkage covariance estimate.

    Parameters
    ----------
    returns : sequence of sequence of float
        ``n`` observations of ``p`` asset returns.

    Returns
    -------
    (sigma_hat, delta) : (list[list[float]], float)
        The shrunk ``p x p`` covariance matrix and the shrinkage intensity
        ``delta`` in ``[0, 1]``. ``delta`` rises toward 1 as the sample estimate
        gets noisier (small ``n``) and falls toward 0 as it gets reliable.
    """
    n, p, means = _demean(returns)
    S = sample_covariance(returns)
    F, r_bar = constant_correlation_target(S)
    std = [math.sqrt(S[i][i]) for i in range(p)]

    dev = [[returns[t][i] - means[i] for i in range(p)] for t in range(n)]

    # pi: sum over (i, j) of Var of the sample covariance entry S_ij.
    pi_mat = [[0.0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            acc = 0.0
            for t in range(n):
                acc += (dev[t][i] * dev[t][j] - S[i][j]) ** 2
            pi_mat[i][j] = acc / n
    pi_hat = sum(sum(row) for row in pi_mat)

    # rho: diagonal terms exactly, off-diagonal via Ledoit-Wolf's asymptotic term.
    rho_hat = sum(pi_mat[i][i] for i in range(p))
    for i in range(p):
        for j in range(p):
            if i == j or std[i] == 0.0 or std[j] == 0.0:
                continue
            # Asymptotic covariance between S_ij and the target entry F_ij.
            theta_ii_ij = 0.0
            theta_jj_ij = 0.0
            for t in range(n):
                d_ii = dev[t][i] * dev[t][i] - S[i][i]
                d_jj = dev[t][j] * dev[t][j] - S[j][j]
                d_ij = dev[t][i] * dev[t][j] - S[i][j]
                theta_ii_ij += d_ii * d_ij
                theta_jj_ij += d_jj * d_ij
            theta_ii_ij /= n
            theta_jj_ij /= n
            rho_hat += 0.5 * r_bar * (std[j] / std[i] * theta_ii_ij
                                      + std[i] / std[j] * theta_jj_ij)

    # gamma: squared Frobenius distance of the target from the sample covariance.
    gamma_hat = 0.0
    for i in range(p):
        for j in range(p):
            gamma_hat += (F[i][j] - S[i][j]) ** 2

    if gamma_hat <= 0.0:
        delta = 0.0
    else:
        kappa = (pi_hat - rho_hat) / gamma_hat
        delta = max(0.0, min(1.0, kappa / n))

    sigma_hat = [[delta * F[i][j] + (1.0 - delta) * S[i][j] for j in range(p)]
                 for i in range(p)]
    return sigma_hat, delta
