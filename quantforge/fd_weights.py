"""Fornberg finite-difference weights on an arbitrary grid.

Given sample locations ``x_0..x_n`` and an evaluation point ``x0``, Fornberg's algorithm
computes the weights ``c_{k,i}`` such that the ``k``-th derivative of a function at ``x0`` is
approximated by ``sum_i c_{k,i} f(x_i)`` -- exact for polynomials up to degree ``n``. The
grid need not be uniform, so this generates any stencil (central, one-sided, staggered) to
any order in one pass. Pure standard library.
"""


def fd_weights(x0, grid, max_deriv):
    """Fornberg weights for derivatives ``0..max_deriv`` at ``x0`` on ``grid``.

    Returns a list ``c`` of length ``max_deriv + 1`` where ``c[k]`` is the weight list (same
    length as ``grid``) for the ``k``-th derivative: ``f^(k)(x0) ~ sum_i c[k][i] f(grid[i])``.
    Requires ``len(grid) > max_deriv`` and distinct grid points.
    """
    n = len(grid) - 1
    if max_deriv < 0:
        raise ValueError("max_deriv must be non-negative")
    if n < max_deriv:
        raise ValueError("need more than max_deriv grid points")
    if len(set(grid)) != len(grid):
        raise ValueError("grid points must be distinct")
    m = max_deriv
    # c[i][j][k] in Fornberg's notation flattened; build weight table
    # weights[k][i]
    c1 = 1.0
    c4 = grid[0] - x0
    # w[i][k]
    w = [[0.0] * (m + 1) for _ in range(n + 1)]
    w[0][0] = 1.0
    for i in range(1, n + 1):
        mn = min(i, m)
        c2 = 1.0
        c5 = c4
        c4 = grid[i] - x0
        for j in range(0, i):
            c3 = grid[i] - grid[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    w[i][k] = c1 * (k * w[i - 1][k - 1] - c5 * w[i - 1][k]) / c2
                w[i][0] = -c1 * c5 * w[i - 1][0] / c2
            for k in range(mn, 0, -1):
                w[j][k] = (c4 * w[j][k] - k * w[j][k - 1]) / c3
            w[j][0] = c4 * w[j][0] / c3
        c1 = c2
    # transpose to weights[k][i]
    return [[w[i][k] for i in range(n + 1)] for k in range(m + 1)]
