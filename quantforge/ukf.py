"""Unscented Kalman filter (UKF) -- derivative-free nonlinear state estimation.

The extended Kalman filter (:mod:`quantforge.ekf`) linearizes ``f`` and ``h`` with their
Jacobians, which can be inaccurate when the nonlinearity is strong. The unscented filter takes a
different route: it deterministically samples ``2n + 1`` *sigma points* around the current mean
(van der Merwe's scaled set), pushes each through the true nonlinear ``f``/``h``, and recovers
the transformed mean and covariance from the propagated points. This captures the mean to
second order and needs no derivatives at all -- ``f`` and ``h`` are ordinary Python functions of
plain floats.

The sigma spread is controlled by ``alpha`` (spread), ``beta`` (prior knowledge; ``2`` is
optimal for Gaussians) and ``kappa`` (secondary scaling). Sigma points use the matrix square
root from :func:`quantforge.linalg.cholesky`. Pure standard library.
"""

import math

from .linalg import cholesky


def _matvec(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def _outer_scaled(a, b, w):
    return [[w * a[i] * b[j] for j in range(len(b))] for i in range(len(a))]


def _add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _inv(A):
    from .lu import lu_solve
    n = len(A)
    cols = [lu_solve(A, [1.0 if i == j else 0.0 for i in range(n)]) for j in range(n)]
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def _sigma_points(x, P, lam):
    n = len(x)
    # matrix square root of (n + lam) P via Cholesky (lower-triangular)
    scaled = [[(n + lam) * P[i][j] for j in range(n)] for i in range(n)]
    L = cholesky(scaled)
    pts = [list(x)]
    for j in range(n):
        col = [L[i][j] for i in range(n)]     # j-th column of the lower factor
        pts.append([x[i] + col[i] for i in range(n)])
        pts.append([x[i] - col[i] for i in range(n)])
    return pts


def unscented_kalman_filter(observations, f, h, Q, R, x0, P0,
                            alpha=1e-3, beta=2.0, kappa=0.0):
    """Unscented Kalman filter over ``observations`` for nonlinear ``f`` and ``h``.

    ``f`` maps a length-``n`` state list to a length-``n`` list (transition); ``h`` maps the
    state to a length-``m`` list (measurement). Both take and return plain floats -- no autodiff
    or Jacobians. ``Q`` (n x n), ``R`` (m x m) covariances; ``x0`` (n), ``P0`` (n x n) initial
    mean/covariance. ``alpha``/``beta``/``kappa`` are the van der Merwe scaling parameters.
    Returns a dict with ``filtered_means`` and ``filtered_covariances``.
    """
    n = len(x0)
    lam = alpha * alpha * (n + kappa) - n
    # weights (van der Merwe)
    wm = [lam / (n + lam)] + [1.0 / (2.0 * (n + lam))] * (2 * n)
    wc = [lam / (n + lam) + (1.0 - alpha * alpha + beta)] + [1.0 / (2.0 * (n + lam))] * (2 * n)

    x = list(x0)
    P = [row[:] for row in P0]
    means, covs = [], []

    for z in observations:
        # ---- predict ----
        sig = _sigma_points(x, P, lam)
        fsig = [f(s) for s in sig]
        xp = [sum(wm[k] * fsig[k][i] for k in range(len(sig))) for i in range(n)]
        Pp = [[Q[i][j] for j in range(n)] for i in range(n)]
        for k in range(len(sig)):
            d = [fsig[k][i] - xp[i] for i in range(n)]
            Pp = _add(Pp, _outer_scaled(d, d, wc[k]))

        # ---- update ----
        sig2 = _sigma_points(xp, Pp, lam)      # re-draw around predicted mean/cov
        hsig = [h(s) for s in sig2]
        m = len(hsig[0])
        zp = [sum(wm[k] * hsig[k][i] for k in range(len(sig2))) for i in range(m)]
        S = [[R[i][j] for j in range(m)] for i in range(m)]
        Cxz = [[0.0] * m for _ in range(n)]
        for k in range(len(sig2)):
            dz = [hsig[k][i] - zp[i] for i in range(m)]
            dx = [sig2[k][i] - xp[i] for i in range(n)]
            S = _add(S, _outer_scaled(dz, dz, wc[k]))
            for i in range(n):
                for j in range(m):
                    Cxz[i][j] += wc[k] * dx[i] * dz[j]
        Sinv = _inv(S)
        # Kalman gain K = Cxz S^{-1}
        K = [[sum(Cxz[i][p] * Sinv[p][j] for p in range(m)) for j in range(m)] for i in range(n)]
        innov = [z[i] - zp[i] for i in range(m)]
        x = [xp[i] + sum(K[i][j] * innov[j] for j in range(m)) for i in range(n)]
        # P = Pp - K S K'
        KS = [[sum(K[i][p] * S[p][j] for p in range(m)) for j in range(m)] for i in range(n)]
        KSKt = [[sum(KS[i][p] * K[j][p] for p in range(m)) for j in range(n)] for i in range(n)]
        P = [[Pp[i][j] - KSKt[i][j] for j in range(n)] for i in range(n)]

        means.append(x[:])
        covs.append([r[:] for r in P])

    return {"filtered_means": means, "filtered_covariances": covs}
