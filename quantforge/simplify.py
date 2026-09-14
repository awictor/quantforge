"""Polyline simplification by the Ramer-Douglas-Peucker algorithm.

Reduces a polyline to fewer points while keeping its shape within a tolerance: the
classic map-generalization and GPS-track-thinning method. It keeps the endpoints, finds
the vertex farthest from the chord between them, and recurses on both halves if that
distance exceeds ``epsilon`` -- otherwise it drops every interior point. Pure standard
library.
"""

import math


def _perp_distance(p, a, b):
    """Perpendicular distance from ``p`` to the segment/line through ``a`` and ``b``."""
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    denom = math.hypot(dx, dy)
    if denom == 0.0:
        return math.hypot(p[0] - ax, p[1] - ay)     # degenerate: distance to the point
    return abs(dx * (ay - p[1]) - dy * (ax - p[0])) / denom


def douglas_peucker(points, epsilon):
    """Simplify a polyline to a subset of its points within perpendicular tolerance ``epsilon``.

    Returns a new list containing the retained points (always including the first and
    last). Larger ``epsilon`` keeps fewer points; ``epsilon = 0`` keeps every point that
    is not exactly collinear. The result is a subsequence of the input in order.
    """
    if epsilon < 0:
        raise ValueError("epsilon must be non-negative")
    n = len(points)
    if n <= 2:
        return list(points)
    # Find the point farthest from the chord (first, last).
    first, last = points[0], points[-1]
    dmax = 0.0
    index = 0
    for i in range(1, n - 1):
        d = _perp_distance(points[i], first, last)
        if d > dmax:
            dmax = d
            index = i
    if dmax > epsilon:
        left = douglas_peucker(points[:index + 1], epsilon)
        right = douglas_peucker(points[index:], epsilon)
        return left[:-1] + right          # drop the duplicated join point
    return [first, last]
