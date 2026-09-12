"""Matrix utilities for risk work: Cholesky, correlation repair, correlated draws.

The Cholesky factorization ``A = L L^T`` of a positive-definite matrix, a check
for positive-definiteness, the transform turning independent standard normals
into correlated ones, and Higham's nearest-positive-semidefinite correlation
matrix (for repairing an indefinite estimated correlation). Pure standard library.
"""

import math


def cholesky(matrix):
    """Lower-triangular Cholesky factor ``L`` with ``L L^T = matrix``.

    Requires a symmetric positive-definite input; raises if a non-positive pivot
    is encountered (matrix not PD).
    """
    n = len(matrix)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                d = matrix[i][i] - s
                if d <= 0.0:
                    raise ValueError("matrix is not positive definite")
                L[i][j] = math.sqrt(d)
            else:
                L[i][j] = (matrix[i][j] - s) / L[j][j]
    return L


def is_positive_definite(matrix):
    """True if the symmetric matrix is positive definite (Cholesky succeeds)."""
    try:
        cholesky(matrix)
        return True
    except ValueError:
        return False


def correlated_normals(independent, correlation):
    """Turn independent standard normals into correlated ones via Cholesky.

    Returns ``L z`` where ``L`` is the Cholesky factor of ``correlation`` -- the
    resulting vector has the target correlation structure. ``independent`` is a
    vector of IID standard normal draws.
    """
    L = cholesky(correlation)
    n = len(L)
    if len(independent) != n:
        raise ValueError("independent vector length must match correlation size")
    return [sum(L[i][k] * independent[k] for k in range(i + 1)) for i in range(n)]


def basket_option_mc(spots, weights, strike, t, r, sigmas, correlation,
                     option_type="call", n_paths=50000, seed=987654321):
    """Monte Carlo basket option on ``sum_i w_i S_i`` for any number of assets.

    Simulates correlated terminal asset prices under geometric Brownian motion
    (correlated normals via :func:`correlated_normals` / Cholesky) and averages the
    discounted basket payoff. A general n-asset reference for the two-asset
    analytic :func:`quantforge.basket_option`. Deterministic per seed.
    """
    n = len(spots)
    if not (len(weights) == len(sigmas) == n):
        raise ValueError("spots, weights, sigmas must have equal length")
    if len(correlation) != n:
        raise ValueError("correlation must be n x n")
    if n_paths < 1:
        raise ValueError("n_paths must be positive")
    is_call = str(option_type).lower() in ("call", "c")
    L = cholesky(correlation)
    disc = math.exp(-r * t)
    sq = math.sqrt(t)
    drift = [(r - 0.5 * s * s) * t for s in sigmas]
    state = seed & 0xFFFFFFFF

    def _normal():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        u1 = (state + 0.5) / 0x80000000
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        u2 = (state + 0.5) / 0x80000000
        return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)

    total = 0.0
    for _ in range(n_paths):
        z = [_normal() for _ in range(n)]
        corr_z = [sum(L[i][k] * z[k] for k in range(i + 1)) for i in range(n)]
        basket = sum(weights[i] * spots[i] * math.exp(drift[i] + sigmas[i] * sq * corr_z[i])
                     for i in range(n))
        payoff = max(basket - strike, 0.0) if is_call else max(strike - basket, 0.0)
        total += payoff
    return disc * total / n_paths


def _sym_eigen(matrix, tol=1e-12, max_sweeps=100):
    from .pca import jacobi_eigen
    return jacobi_eigen(matrix, tol, max_sweeps)


def nearest_correlation(matrix, max_iter=100, tol=1e-10):
    """Nearest positive-semidefinite correlation matrix (Higham alternating projection).

    Repairs an indefinite estimated correlation matrix to the closest valid one:
    alternately projects onto the PSD cone (clip negative eigenvalues to zero) and
    onto the unit-diagonal set, iterating to convergence. Returns a symmetric PSD
    matrix with unit diagonal; a matrix that is already a valid correlation is
    returned essentially unchanged.
    """
    n = len(matrix)
    Y = [row[:] for row in matrix]
    dS = [[0.0] * n for _ in range(n)]
    for _ in range(max_iter):
        R = [[Y[i][j] - dS[i][j] for j in range(n)] for i in range(n)]
        # Project R onto the PSD cone via eigenvalue clipping.
        vals, vecs = _sym_eigen(R)
        X = [[0.0] * n for _ in range(n)]
        for k in range(n):
            lam = max(vals[k], 0.0)
            if lam == 0.0:
                continue
            vk = vecs[k]
            for i in range(n):
                for j in range(n):
                    X[i][j] += lam * vk[i] * vk[j]
        dS = [[X[i][j] - R[i][j] for j in range(n)] for i in range(n)]
        # Project onto unit-diagonal set.
        newY = [[X[i][j] for j in range(n)] for i in range(n)]
        for i in range(n):
            newY[i][i] = 1.0
        diff = max(abs(newY[i][j] - Y[i][j]) for i in range(n) for j in range(n))
        Y = newY
        if diff < tol:
            break
    return Y
