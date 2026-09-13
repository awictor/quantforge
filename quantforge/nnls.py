"""Non-negative least squares (Lawson-Hanson active-set algorithm).

Solve ``min ||A x - b||^2`` subject to ``x >= 0``. The non-negativity constraint is
natural whenever coefficients are physically weights or proportions -- unmixing spectra,
fitting factor exposures, portfolio-style attribution -- and it often yields a sparse
solution (many coefficients exactly zero) without an explicit penalty. The classic
Lawson-Hanson (1974) active-set method solves it exactly in finitely many steps:
maintain a passive set of free variables, solve the unconstrained least squares on it,
and swap variables in or out to enforce ``x >= 0``. Pure standard library.
"""

from .portopt import _invert


def _lstsq_on_set(A, b, cols):
    """Unconstrained least squares of ``b`` on the columns ``cols`` of ``A``."""
    m = len(A)
    k = len(cols)
    # Normal equations on the sub-matrix.
    ata = [[sum(A[r][cols[i]] * A[r][cols[j]] for r in range(m)) for j in range(k)]
           for i in range(k)]
    atb = [sum(A[r][cols[i]] * b[r] for r in range(m)) for i in range(k)]
    inv = _invert(ata)
    return [sum(inv[i][j] * atb[j] for j in range(k)) for i in range(k)]


def nnls(A, b, max_iter=None, tol=1e-10):
    """Non-negative least squares: ``min ||A x - b||^2`` with ``x >= 0``.

    ``A`` is an ``m x n`` matrix (list of rows), ``b`` a length-``m`` vector. Returns a
    dict with the non-negative solution ``x`` (length ``n``), the ``residual_norm``
    ``||A x - b||``, and ``n_iter``. Uses the Lawson-Hanson active-set method.
    """
    m = len(A)
    if m != len(b):
        raise ValueError("A and b must have matching row counts")
    if m == 0:
        raise ValueError("need at least one row")
    n = len(A[0])
    if max_iter is None:
        max_iter = 3 * n

    x = [0.0] * n
    passive = []                      # indices with x free (currently > 0)
    active = list(range(n))           # indices constrained to 0

    def gradient(x):
        # w = A'(b - A x): the negative gradient of the objective.
        resid = [b[r] - sum(A[r][j] * x[j] for j in range(n)) for r in range(m)]
        return [sum(A[r][j] * resid[r] for r in range(m)) for j in range(n)]

    it = 0
    while active and it < max_iter:
        w = gradient(x)
        # Pick the active index with the largest positive gradient (most promising).
        j_star = max(active, key=lambda j: w[j])
        if w[j_star] <= tol:
            break                     # KKT satisfied: no active var improves the fit
        active.remove(j_star)
        passive.append(j_star)

        # Inner loop: solve LS on passive set, drop any non-positive entries.
        while True:
            it += 1
            if it > max_iter:
                break
            passive.sort()
            z_vals = _lstsq_on_set(A, b, passive)
            z = dict(zip(passive, z_vals))
            if all(z[j] > tol for j in passive):
                for j in passive:
                    x[j] = z[j]
                break
            # Step toward z only as far as non-negativity allows.
            alpha = min(x[j] / (x[j] - z[j]) for j in passive
                        if z[j] <= tol and x[j] - z[j] != 0.0)
            x = [x[j] + alpha * (z.get(j, 0.0) - x[j]) for j in range(n)]
            # Move variables that hit zero back to the active set.
            newly = [j for j in passive if abs(x[j]) <= tol]
            for j in newly:
                passive.remove(j)
                active.append(j)
                x[j] = 0.0

    resid = [b[r] - sum(A[r][j] * x[j] for j in range(n)) for r in range(m)]
    rn = sum(rr * rr for rr in resid) ** 0.5
    return {"x": x, "residual_norm": rn, "n_iter": it}
