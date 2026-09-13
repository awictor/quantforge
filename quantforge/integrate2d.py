"""Two-dimensional numerical integration over a rectangle.

Extends the 1-D quadrature to double integrals ``int_a^b int_c^d f(x, y) dy dx``:

  * ``integrate2d_gauss`` -- a tensor-product Gauss-Legendre rule (order n in each
    axis), exact for polynomials up to degree ``2n-1`` per variable and very accurate
    for smooth integrands with only ``n^2`` evaluations.
  * ``integrate2d_simpson`` -- composite Simpson on a grid, for less-smooth integrands.

Both integrate over an axis-aligned rectangle. Pure standard library; reuses the
library's Gauss-Legendre nodes.
"""

from .quadrature import _GL, simpson


def integrate2d_gauss(f, a, b, c, d, n=5):
    """Tensor Gauss-Legendre integral of ``f(x, y)`` over ``[a, b] x [c, d]``.

    ``n`` in {2, 3, 4, 5} is the order per axis; exact for polynomials up to degree
    ``2n - 1`` in each variable. Uses ``n^2`` function evaluations.
    """
    if n not in _GL:
        raise ValueError("n must be one of 2, 3, 4, 5")
    nodes, weights = _GL[n]
    hx = 0.5 * (b - a)
    mx = 0.5 * (a + b)
    hy = 0.5 * (d - c)
    my = 0.5 * (c + d)
    total = 0.0
    for i in range(n):
        xi = mx + hx * nodes[i]
        for j in range(n):
            yj = my + hy * nodes[j]
            total += weights[i] * weights[j] * f(xi, yj)
    return hx * hy * total


def integrate2d_simpson(f, a, b, c, d, nx=50, ny=50):
    """Composite-Simpson double integral of ``f(x, y)`` over ``[a, b] x [c, d]``.

    Integrates in ``y`` for each ``x`` then in ``x`` (Fubini), each by composite
    Simpson with ``ny`` / ``nx`` panels. Slower but robust for integrands the smooth
    Gauss rule would miss.
    """
    def inner(x):
        return simpson(lambda y: f(x, y), c, d, ny)
    return simpson(inner, a, b, nx)
