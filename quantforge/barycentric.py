"""Barycentric Lagrange interpolation: stable, reusable polynomial interpolation.

Neville's algorithm re-does ``O(n^2)`` work for every evaluation point. The
barycentric form precomputes weights once, then evaluates the interpolant at any point
in ``O(n)`` with the formula

    p(x) = [ sum_j w_j y_j / (x - x_j) ] / [ sum_j w_j / (x - x_j) ],

which is numerically stable (unlike the naive Lagrange or monomial forms) and lets you
add points or change the sampled values cheaply. On Chebyshev nodes the weights have a
closed form and interpolation converges spectrally, avoiding the Runge blow-up of
equally-spaced nodes. Pure standard library.
"""

import math


def barycentric_weights(xs):
    """Barycentric weights ``w_j = 1 / prod_{k != j} (x_j - x_k)`` for nodes ``xs``.

    ``O(n^2)`` one-time cost; feed the result to :func:`barycentric_eval`. Raises on
    duplicate nodes.
    """
    n = len(xs)
    if n == 0:
        raise ValueError("need at least one node")
    w = [1.0] * n
    for j in range(n):
        for k in range(n):
            if k != j:
                diff = xs[j] - xs[k]
                if diff == 0.0:
                    raise ValueError("duplicate nodes")
                w[j] /= diff
    return w


def barycentric_eval(xs, ys, weights, x):
    """Evaluate the barycentric interpolant at ``x`` in ``O(n)``.

    ``weights`` come from :func:`barycentric_weights`. Handles the case where ``x``
    coincides with a node exactly (returns that node's value).
    """
    n = len(xs)
    if not (n == len(ys) == len(weights)):
        raise ValueError("xs, ys, weights must have equal length")
    num = 0.0
    den = 0.0
    for j in range(n):
        d = x - xs[j]
        if d == 0.0:
            return ys[j]
        t = weights[j] / d
        num += t * ys[j]
        den += t
    if den == 0.0:
        raise ValueError("barycentric denominator vanished")
    return num / den


def chebyshev_nodes(a, b, n):
    """The ``n`` Chebyshev points of the second kind (Chebyshev-Lobatto) on ``[a, b]``.

    ``x_i = (a+b)/2 + (b-a)/2 * cos(i pi / (n-1))`` -- the nodes that make polynomial
    interpolation converge spectrally for smooth functions and avoid Runge
    oscillation. Returns them in increasing order.
    """
    if n < 2:
        raise ValueError("need at least 2 nodes")
    mid = 0.5 * (a + b)
    half = 0.5 * (b - a)
    pts = [mid + half * math.cos(i * math.pi / (n - 1)) for i in range(n)]
    pts.sort()
    return pts


def chebyshev_barycentric_weights(n):
    """Closed-form barycentric weights for the ``n`` Chebyshev-Lobatto nodes.

    ``w_i = (-1)^i delta_i`` with the endpoints halved (``delta = 1/2`` at ``i = 0,
    n-1`` and ``1`` inside). These pair with :func:`chebyshev_nodes` and are far more
    stable than recomputing the general product. Returned in increasing-``x`` order to
    match :func:`chebyshev_nodes`.
    """
    if n < 2:
        raise ValueError("need at least 2 nodes")
    # cos(i pi/(n-1)) is decreasing in i, so increasing-x order reverses the index.
    w = []
    for i in range(n):
        wi = 1.0 if 0 < i < n - 1 else 0.5
        w.append(((-1) ** i) * wi)
    w.reverse()
    return w
