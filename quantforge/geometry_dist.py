"""Point-to-line and point-to-segment distances in the plane.

The perpendicular distance from a point to an infinite line, the (clamped) distance to a
finite segment, the closest point on a segment, and the nearest distance from a point to a
polyline. The workhorse for hit-testing, snapping, and route proximity. Points are
``(x, y)`` tuples. Pure standard library.
"""

import math


def point_to_line_distance(p, a, b):
    """Perpendicular distance from point ``p`` to the *infinite* line through ``a`` and ``b``.

    Raises if ``a == b`` (no line defined).
    """
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    denom = math.hypot(dx, dy)
    if denom == 0.0:
        raise ValueError("a and b must be distinct")
    # |cross product| / |direction|.
    return abs(dx * (ay - p[1]) - dy * (ax - p[0])) / denom


def closest_point_on_segment(p, a, b):
    """Closest point on the finite segment ``a-b`` to ``p`` (clamped to the endpoints)."""
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    len2 = dx * dx + dy * dy
    if len2 == 0.0:
        return (ax, ay)                 # degenerate segment = a point
    t = ((p[0] - ax) * dx + (p[1] - ay) * dy) / len2
    t = max(0.0, min(1.0, t))
    return (ax + t * dx, ay + t * dy)


def point_segment_distance(p, a, b):
    """Distance from ``p`` to the nearest point of the finite segment ``a-b``.

    Equals the perpendicular distance when the foot of the perpendicular lands on the
    segment, otherwise the distance to the nearer endpoint.
    """
    cx, cy = closest_point_on_segment(p, a, b)
    return math.hypot(p[0] - cx, p[1] - cy)


def point_polyline_distance(p, polyline):
    """Minimum distance from ``p`` to a polyline (a list of ``>= 2`` vertices).

    Takes the smallest :func:`point_segment_distance` over consecutive segments.
    """
    if len(polyline) < 2:
        raise ValueError("polyline needs at least two vertices")
    best = float("inf")
    for i in range(len(polyline) - 1):
        d = point_segment_distance(p, polyline[i], polyline[i + 1])
        if d < best:
            best = d
    return best
