"""Two-dimensional Poisson equation solver (successive over-relaxation).

Solves the elliptic PDE ``u_xx + u_yy = f(x, y)`` on a rectangular grid with Dirichlet
boundary values, by successive over-relaxation (SOR) -- Gauss-Seidel sweeps accelerated
by an over-relaxation factor ``omega``. With ``f = 0`` this is the Laplace equation, and
its solution is harmonic (the boundary values diffuse smoothly into the interior with no
interior extrema). Iterates to a residual tolerance. Pure standard library.
"""


def poisson2d(f_grid, boundary, dx, dy, omega=1.5, tol=1e-8, max_iter=10000):
    """Solve ``u_xx + u_yy = f`` on a grid by SOR with Dirichlet boundaries.

    ``f_grid`` is the source term ``f`` as a 2-D list (rows = y, cols = x);
    ``boundary`` a same-shaped grid whose *edge* values set the fixed boundary (interior
    entries are the initial guess). ``dx``/``dy`` are the grid spacings, ``omega`` the
    SOR factor in ``(0, 2)`` (1 = Gauss-Seidel). Returns ``(u, n_iter)`` -- the solution
    grid and the sweeps taken. Laplace is the ``f_grid`` all-zero case.
    """
    ny = len(f_grid)
    nx = len(f_grid[0])
    if ny < 3 or nx < 3:
        raise ValueError("need at least a 3x3 grid")
    if not (0.0 < omega < 2.0):
        raise ValueError("omega must be in (0, 2)")
    u = [list(row) for row in boundary]
    hx2 = dx * dx
    hy2 = dy * dy
    denom = 2.0 * (hx2 + hy2)

    n_iter = 0
    for it in range(1, max_iter + 1):
        n_iter = it
        max_resid = 0.0
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                # Gauss-Seidel update of the 5-point Laplacian = f.
                new = (hy2 * (u[j][i - 1] + u[j][i + 1])
                       + hx2 * (u[j - 1][i] + u[j + 1][i])
                       - hx2 * hy2 * f_grid[j][i]) / denom
                delta = new - u[j][i]
                u[j][i] += omega * delta
                if abs(delta) > max_resid:
                    max_resid = abs(delta)
        if max_resid < tol:
            break
    return u, n_iter
