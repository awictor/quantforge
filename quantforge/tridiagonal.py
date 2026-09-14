"""Tridiagonal and cyclic-tridiagonal linear solvers.

A tridiagonal system -- nonzero entries only on the main diagonal and its two neighbours --
arises in cubic splines, implicit PDE steps, and 1-D diffusion. The Thomas algorithm solves
it in ``O(n)`` instead of the ``O(n^3)`` of a dense solve. The *cyclic* variant, where the
first and last equations also couple (periodic boundary), is handled by the Sherman-Morrison
correction over two Thomas solves. Pure standard library.
"""


def solve_tridiagonal(lower, diag, upper, rhs):
    """Solve a tridiagonal system by the Thomas algorithm; returns the solution list.

    ``diag`` is the main diagonal (length ``n``); ``lower[i]`` is the sub-diagonal entry in
    row ``i`` (``lower[0]`` unused); ``upper[i]`` the super-diagonal entry in row ``i``
    (``upper[n-1]`` unused). All four lists have length ``n``. Raises on a zero pivot (the
    matrix is singular or not diagonally dominant enough for plain elimination).
    """
    n = len(diag)
    if not (len(lower) == len(upper) == len(rhs) == n):
        raise ValueError("all inputs must have the same length")
    if n == 0:
        return []
    c = [0.0] * n      # modified super-diagonal
    d = [0.0] * n      # modified rhs
    if diag[0] == 0:
        raise ValueError("zero pivot at row 0")
    c[0] = upper[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, n):
        m = diag[i] - lower[i] * c[i - 1]
        if m == 0:
            raise ValueError("zero pivot at row %d" % i)
        c[i] = upper[i] / m
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / m
    x = [0.0] * n
    x[n - 1] = d[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


def solve_cyclic_tridiagonal(lower, diag, upper, rhs):
    """Solve a *cyclic* tridiagonal system (periodic boundary) via Sherman-Morrison.

    Same layout as :func:`solve_tridiagonal`, but ``lower[0]`` is the corner entry coupling
    row 0 to column ``n-1``, and ``upper[n-1]`` couples row ``n-1`` to column 0. Uses two
    Thomas solves plus a rank-1 correction. Needs ``n >= 3``.
    """
    n = len(diag)
    if not (len(lower) == len(upper) == len(rhs) == n):
        raise ValueError("all inputs must have the same length")
    if n < 3:
        raise ValueError("cyclic solve needs n >= 3")
    alpha = lower[0]        # corner (0, n-1)
    beta = upper[n - 1]     # corner (n-1, 0)
    # choose gamma to avoid a zero on the modified first pivot
    gamma = -diag[0] if diag[0] != 0 else 1.0
    b = list(diag)
    b[0] = diag[0] - gamma
    b[n - 1] = diag[n - 1] - alpha * beta / gamma
    # solve A' x = rhs
    x = solve_tridiagonal(lower, b, upper, rhs)
    # solve A' z = u, where u = gamma e_0 + beta e_{n-1}
    u = [0.0] * n
    u[0] = gamma
    u[n - 1] = beta
    z = solve_tridiagonal(lower, b, upper, u)
    fact = (x[0] + alpha * x[n - 1] / gamma) / (1.0 + z[0] + alpha * z[n - 1] / gamma)
    return [x[i] - fact * z[i] for i in range(n)]
