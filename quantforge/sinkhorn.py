"""Entropic optimal transport via the Sinkhorn-Knopp algorithm.

Optimal transport measures the cost of morphing one distribution into another: given a cost
matrix ``C`` (``C[i][j]`` = cost of moving a unit of mass from source ``i`` to target ``j``) and
marginals ``a``, ``b``, find the transport plan ``P >= 0`` with row sums ``a`` and column sums
``b`` minimizing ``sum P[i][j] C[i][j]``. The exact problem is a linear program; adding an
entropy penalty ``-eps H(P)`` makes it strongly convex and solvable by the fast matrix-scaling
iteration of Sinkhorn-Knopp -- alternately rescaling rows to match ``a`` and columns to match
``b``. Smaller ``eps`` approaches the true (unregularized) transport cost.

This underlies distribution comparison (a multivariate generalization of the 1-D
:func:`quantforge.wasserstein.wasserstein_distance`), color transfer, and domain adaptation.
Pure standard library.
"""

import math


def _lse(vals):
    # numerically stable log-sum-exp of a list
    mx = max(vals)
    if mx == float("-inf"):
        return float("-inf")
    return mx + math.log(sum(math.exp(v - mx) for v in vals))


def sinkhorn(a, b, C, eps=0.05, max_iter=2000, tol=1e-9):
    """Entropic-regularized optimal transport plan between marginals ``a`` and ``b``.

    ``a`` (length ``n``) and ``b`` (length ``m``) are non-negative weight vectors with equal
    total mass; ``C`` is the ``n x m`` cost matrix; ``eps`` the entropic regularization strength.
    Returns a dict with ``plan`` (the ``n x m`` transport matrix), ``cost`` (``sum P*C``) and
    ``n_iter``. As ``eps -> 0`` the cost approaches the exact optimal-transport cost.

    Runs in the log-domain (stabilized potentials ``f, g``) so it stays accurate for small
    ``eps`` where the raw ``exp(-C/eps)`` kernel would underflow to zero.
    """
    n, m = len(a), len(b)
    sa, sb = sum(a), sum(b)
    if abs(sa - sb) > 1e-6 * max(1.0, sa):
        raise ValueError("marginals a and b must have equal total mass")
    log_a = [math.log(ai) if ai > 0 else float("-inf") for ai in a]
    log_b = [math.log(bj) if bj > 0 else float("-inf") for bj in b]
    # dual potentials; plan P[i][j] = exp((f[i] + g[j] - C[i][j]) / eps)
    f = [0.0] * n
    g = [0.0] * m
    n_iter = 0
    for it in range(max_iter):
        n_iter = it + 1
        f_prev = f[:]
        # f_i = eps*(log a_i - LSE_j((g_j - C_ij)/eps))
        for i in range(n):
            f[i] = eps * (log_a[i] - _lse([(g[j] - C[i][j]) / eps for j in range(m)]))
        # g_j = eps*(log b_j - LSE_i((f_i - C_ij)/eps))
        for j in range(m):
            g[j] = eps * (log_b[j] - _lse([(f[i] - C[i][j]) / eps for i in range(n)]))
        err = max(abs(f[i] - f_prev[i]) for i in range(n))
        if err < tol:
            break
    plan = [[math.exp((f[i] + g[j] - C[i][j]) / eps) for j in range(m)] for i in range(n)]
    cost = sum(plan[i][j] * C[i][j] for i in range(n) for j in range(m))
    return {"plan": plan, "cost": cost, "n_iter": n_iter}


def cost_matrix(xs, ys, p=2):
    """Ground-cost matrix ``C[i][j] = |xs[i] - ys[j]|^p`` for scalar or vector support points.

    ``xs`` and ``ys`` are lists of points (each a float or an equal-length coordinate list);
    ``p`` is the exponent of the Euclidean distance (``p = 2`` gives squared distance, the
    2-Wasserstein ground cost). Returns the ``len(xs) x len(ys)`` matrix.
    """
    def dist(a, b):
        if isinstance(a, (list, tuple)):
            return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))
        return abs(a - b)

    return [[dist(x, y) ** p for y in ys] for x in xs]
