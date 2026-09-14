"""Bezier curves and Bernstein polynomials (de Casteljau's algorithm).

A Bezier curve of degree ``n`` is the Bernstein-weighted blend of ``n + 1`` control points,
``B(t) = sum_i C(n, i) t^i (1-t)^(n-i) P_i`` for ``t`` in ``[0, 1]``. De Casteljau's
algorithm evaluates it by repeated linear interpolation -- numerically stable and, as a
by-product, splits the curve at ``t`` into two sub-curves (subdivision). This module
evaluates a curve, its derivative (a degree ``n-1`` Bezier of the control differences),
subdivides it, and exposes the Bernstein basis directly. Control points are scalars or
equal-length coordinate tuples. Pure standard library.
"""

from math import comb


def _lerp(p, q, t):
    if isinstance(p, (int, float)):
        return p + (q - p) * t
    return tuple(pi + (qi - pi) * t for pi, qi in zip(p, q))


def bernstein(n, i, t):
    """Bernstein basis polynomial ``b_{i,n}(t) = C(n, i) t^i (1-t)^(n-i)``."""
    if not 0 <= i <= n:
        return 0.0
    return comb(n, i) * (t ** i) * ((1 - t) ** (n - i))


def bezier_point(control, t):
    """Evaluate the Bezier curve with the given ``control`` points at parameter ``t``.

    Uses de Casteljau's repeated-interpolation scheme. Control points may be scalars or
    coordinate tuples (all the same length).
    """
    if not control:
        raise ValueError("need at least one control point")
    pts = list(control)
    while len(pts) > 1:
        pts = [_lerp(pts[i], pts[i + 1], t) for i in range(len(pts) - 1)]
    return pts[0]


def bezier_curve(control, samples):
    """Sample the Bezier curve at ``samples`` equally spaced ``t`` in ``[0, 1]`` (inclusive)."""
    if samples < 2:
        raise ValueError("samples must be at least 2")
    return [bezier_point(control, i / (samples - 1)) for i in range(samples)]


def bezier_derivative_control(control):
    """Control points of the curve's derivative: a degree ``n-1`` Bezier.

    The derivative of a degree-``n`` Bezier is degree ``n-1`` with control points
    ``n (P_{i+1} - P_i)``. Returns those points (empty for a single control point).
    """
    n = len(control) - 1
    if n < 1:
        return []
    out = []
    for i in range(n):
        p, q = control[i], control[i + 1]
        if isinstance(p, (int, float)):
            out.append(n * (q - p))
        else:
            out.append(tuple(n * (qi - pi) for pi, qi in zip(p, q)))
    return out


def bezier_tangent(control, t):
    """The tangent vector (derivative) of the curve at ``t``."""
    d = bezier_derivative_control(control)
    if not d:
        # a single point has zero derivative
        p = control[0]
        return 0.0 if isinstance(p, (int, float)) else tuple(0.0 for _ in p)
    return bezier_point(d, t)


def bezier_subdivide(control, t):
    """Split the curve at ``t`` into two control-point lists (left, right).

    The de Casteljau intermediate points give both halves exactly: the left curve uses the
    first point of each interpolation level, the right curve the last.
    """
    pts = list(control)
    left = [pts[0]]
    right = [pts[-1]]
    while len(pts) > 1:
        pts = [_lerp(pts[i], pts[i + 1], t) for i in range(len(pts) - 1)]
        left.append(pts[0])
        right.append(pts[-1])
    right.reverse()
    return left, right
