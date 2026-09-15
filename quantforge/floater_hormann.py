"""Floater-Hormann barycentric rational interpolation (pole-free, Runge-stable).

Polynomial interpolation on equispaced nodes suffers the Runge phenomenon: the error blows up
near the endpoints as the degree grows. Barycentric *rational* interpolation with the
Floater-Hormann weights avoids it. For a blending parameter ``d`` (0 <= d <= n), the interpolant
is a blend of local degree-``d`` polynomials; it interpolates every node exactly, has **no poles**
anywhere on the real line, and approximates smooth functions on equispaced nodes with
``O(h^{d+1})`` accuracy -- no Runge blow-up.

This computes the Floater-Hormann weights and evaluates the interpolant via the barycentric
formula (reusing :func:`quantforge.barycentric.barycentric_eval`). Pure standard library.
"""

from .barycentric import barycentric_eval


def floater_hormann_weights(xs, d=3):
    """Floater-Hormann barycentric weights for nodes ``xs`` with blending degree ``d``.

    ``w_k = sum_{i in J_k} (-1)^i / prod_{j != k, j in [i, i+d]} (x_k - x_j)`` where ``J_k`` runs
    over the local index windows containing ``k``. ``d = 0`` gives the Berrut (pole-free) rational
    interpolant; larger ``d`` raises the approximation order. Requires ``0 <= d <= len(xs) - 1``.
    """
    n = len(xs)
    if n == 0:
        raise ValueError("need at least one node")
    if not (0 <= d <= n - 1):
        raise ValueError("require 0 <= d <= len(xs) - 1")
    w = [0.0] * n
    for k in range(n):
        # windows [i, i+d] that contain node k, with window start i in [0, n-1-d]
        i_lo = max(0, k - d)
        i_hi = min(k, n - 1 - d)
        s = 0.0
        for i in range(i_lo, i_hi + 1):
            prod = 1.0
            for j in range(i, i + d + 1):
                if j != k:
                    prod *= (xs[k] - xs[j])
            s += ((-1.0) ** i) / prod
        w[k] = s
    return w


def floater_hormann_interpolate(xs, ys, x, d=3):
    """Evaluate the Floater-Hormann rational interpolant through ``(xs, ys)`` at ``x``.

    Convenience wrapper: builds the weights for blending degree ``d`` and evaluates via the
    barycentric formula. Interpolates the nodes exactly and never has a pole.
    """
    w = floater_hormann_weights(xs, d)
    return barycentric_eval(xs, ys, w, x)
