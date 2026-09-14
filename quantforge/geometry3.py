"""Planar geometry III: diameter, bounding box, minimum enclosing circle.

Extremal measures of a point set: the diameter (farthest pair), the axis-aligned
bounding box, and the smallest circle containing every point (Welzl's expected-linear
algorithm). These answer "how big is this cloud" for spatial indexing, collision bounds,
and clustering. Pure standard library on top of :mod:`quantforge.geometry`.
"""

import math

from .geometry import convex_hull


def bounding_box(points):
    """Axis-aligned bounding box of a point set as ``(min_x, min_y, max_x, max_y)``."""
    if not points:
        raise ValueError("need at least one point")
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return (min(xs), min(ys), max(xs), max(ys))


def polygon_diameter(points):
    """Diameter of a point set: the farthest-apart pair and their distance.

    Returns ``(p, q, distance)``. Restricts the search to the convex-hull vertices (the
    diameter is always realized by two hull points), then does an all-pairs scan over the
    hull -- exact, and cheap once the hull has few vertices.
    """
    if len(points) < 2:
        raise ValueError("need at least two points")
    hull = convex_hull(points)
    if len(hull) < 2:
        # all points coincide
        return (tuple(points[0]), tuple(points[0]), 0.0)
    best = None
    for i in range(len(hull)):
        for j in range(i + 1, len(hull)):
            d = math.hypot(hull[i][0] - hull[j][0], hull[i][1] - hull[j][1])
            if best is None or d > best[2]:
                best = (hull[i], hull[j], d)
    return best


def _circle_from_two(a, b):
    cx = (a[0] + b[0]) / 2.0
    cy = (a[1] + b[1]) / 2.0
    r = math.hypot(a[0] - b[0], a[1] - b[1]) / 2.0
    return (cx, cy, r)


def _circle_from_three(a, b, c):
    ax, ay = a
    bx, by = b
    cx_, cy_ = c
    d = 2.0 * (ax * (by - cy_) + bx * (cy_ - ay) + cx_ * (ay - by))
    if abs(d) < 1e-18:
        return None                      # collinear
    ux = ((ax * ax + ay * ay) * (by - cy_) + (bx * bx + by * by) * (cy_ - ay) +
          (cx_ * cx_ + cy_ * cy_) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx_ - bx) + (bx * bx + by * by) * (ax - cx_) +
          (cx_ * cx_ + cy_ * cy_) * (bx - ax)) / d
    r = math.hypot(ux - ax, uy - ay)
    return (ux, uy, r)


def _in_circle(circle, p, tol=1e-9):
    return math.hypot(p[0] - circle[0], p[1] - circle[1]) <= circle[2] + tol


def min_enclosing_circle(points, seed=1234567):
    """Smallest circle containing all ``points`` as ``(center_x, center_y, radius)``.

    Welzl's algorithm with a deterministic shuffle (an LCG seeded by ``seed``), running in
    expected linear time. Every input point lies inside or on the returned circle, and the
    circle is the unique smallest such -- no smaller circle contains them all.
    """
    if not points:
        raise ValueError("need at least one point")
    pts = [tuple(p) for p in points]
    # Deterministic Fisher-Yates shuffle so the result is reproducible.
    state = seed & 0x7FFFFFFF or 1
    for i in range(len(pts) - 1, 0, -1):
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        j = state % (i + 1)
        pts[i], pts[j] = pts[j], pts[i]

    circle = None
    for i, p in enumerate(pts):
        if circle is None or not _in_circle(circle, p):
            circle = (p[0], p[1], 0.0)
            for j in range(i):
                q = pts[j]
                if not _in_circle(circle, q):
                    circle = _circle_from_two(p, q)
                    for k in range(j):
                        s = pts[k]
                        if not _in_circle(circle, s):
                            c3 = _circle_from_three(p, q, s)
                            if c3 is not None:
                                circle = c3
    return circle
