"""Chebyshev polynomial approximation of a function on an interval.

Approximates ``f`` on ``[a, b]`` by a truncated Chebyshev series
``f(x) ~ sum_k c_k T_k(u)`` with ``u = (2x - a - b)/(b - a)`` in ``[-1, 1]``. For a
smooth function the coefficients decay geometrically, so a modest degree gives
near-machine-precision accuracy -- the basis behind spectral methods and fast
function evaluation. The coefficients are obtained from samples at the Chebyshev
extrema (a discrete cosine transform); evaluation uses the numerically stable
Clenshaw recurrence. Pure standard library.
"""

import math


def chebyshev_fit(f, a, b, degree):
    """Chebyshev coefficients of ``f`` on ``[a, b]`` up to ``degree``.

    Samples ``f`` at the ``degree + 1`` Chebyshev extrema and takes the discrete
    cosine transform. Returns the coefficient list ``[c_0, ..., c_degree]`` (the
    ``c_0`` term is the mean level, used at half weight by :func:`chebyshev_eval`).
    Exact for polynomials of degree ``<= degree``.
    """
    if degree < 1:
        raise ValueError("degree must be at least 1")
    n = degree
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    # Function values at x_j = mid + half*cos(pi j / n).
    fx = [f(mid + half * math.cos(math.pi * j / n)) for j in range(n + 1)]
    coeffs = []
    for k in range(n + 1):
        s = 0.0
        for j in range(n + 1):
            w = 0.5 if (j == 0 or j == n) else 1.0
            s += w * fx[j] * math.cos(math.pi * k * j / n)
        coeffs.append((2.0 / n) * s)
    # On the Chebyshev-Lobatto grid the top coefficient carries a half weight too
    # (mirror of the c_0 half weight); without it a degree-n polynomial aliases.
    coeffs[n] *= 0.5
    return coeffs


def chebyshev_eval(coeffs, a, b, x):
    """Evaluate a Chebyshev series at ``x`` by the Clenshaw recurrence.

    ``coeffs`` are from :func:`chebyshev_fit`; the ``c_0`` term is taken at half
    weight (the standard convention). Stable across ``[a, b]``.
    """
    n = len(coeffs)
    if n == 0:
        raise ValueError("need at least one coefficient")
    u = (2.0 * x - a - b) / (b - a)
    u2 = 2.0 * u
    d = 0.0
    dd = 0.0
    for k in range(n - 1, 0, -1):
        d, dd = u2 * d - dd + coeffs[k], d
    return u * d - dd + 0.5 * coeffs[0]


def chebyshev_derivative(coeffs, a, b):
    """Chebyshev coefficients of the derivative of a series on ``[a, b]``.

    Applies the standard Chebyshev differentiation recurrence and the chain-rule
    factor ``2 / (b - a)``. The returned coefficients evaluate (via
    :func:`chebyshev_eval`) to ``f'`` on the same interval.
    """
    n = len(coeffs)
    if n < 2:
        return [0.0]
    # Derivative coefficients on [-1, 1] (recurrence downward).
    dc = [0.0] * n
    dc[n - 1] = 0.0
    if n >= 2:
        dc[n - 2] = 2.0 * (n - 1) * coeffs[n - 1]
    for k in range(n - 3, -1, -1):
        dc[k] = dc[k + 2] + 2.0 * (k + 1) * coeffs[k + 1]
    scale = 2.0 / (b - a)
    return [scale * c for c in dc[:n - 1]]
