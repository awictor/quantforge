"""Circle geometry: circumcircle, and circle-line / circle-circle intersection.

The circle through three points (circumcircle), where a circle meets an infinite line,
and where two circles meet -- the primitives behind geometric construction, collision
detection, and trilateration. Circles are ``(cx, cy, r)``; points/lines are ``(x, y)``
tuples. Pure standard library.
"""

import math


def circle_from_3points(a, b, c):
    """Circle ``(cx, cy, r)`` through three non-collinear points.

    Raises ``ValueError`` if the points are collinear (no finite circle).
    """
    ax, ay = a
    bx, by = b
    cx_, cy_ = c
    d = 2.0 * (ax * (by - cy_) + bx * (cy_ - ay) + cx_ * (ay - by))
    if abs(d) < 1e-15:
        raise ValueError("points are collinear; no circumcircle")
    ux = ((ax * ax + ay * ay) * (by - cy_) + (bx * bx + by * by) * (cy_ - ay) +
          (cx_ * cx_ + cy_ * cy_) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx_ - bx) + (bx * bx + by * by) * (ax - cx_) +
          (cx_ * cx_ + cy_ * cy_) * (bx - ax)) / d
    r = math.hypot(ux - ax, uy - ay)
    return (ux, uy, r)


def point_in_circle(p, circle):
    """Whether point ``p`` lies inside or on the circle ``(cx, cy, r)``."""
    cx, cy, r = circle
    return math.hypot(p[0] - cx, p[1] - cy) <= r + 1e-12


def circle_line_intersection(circle, a, b):
    """Intersection points of a circle with the *infinite* line through ``a`` and ``b``.

    Returns a list of 0, 1 (tangent), or 2 ``(x, y)`` points. ``a`` and ``b`` must be
    distinct.
    """
    cx, cy, r = circle
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    if dx == 0.0 and dy == 0.0:
        raise ValueError("a and b must be distinct")
    # Project the center onto the line; solve the quadratic in the line parameter.
    fx, fy = ax - cx, ay - cy
    aa = dx * dx + dy * dy
    bb = 2.0 * (fx * dx + fy * dy)
    cc = fx * fx + fy * fy - r * r
    disc = bb * bb - 4.0 * aa * cc
    if disc < 0:
        return []
    if abs(disc) < 1e-12:
        t = -bb / (2.0 * aa)
        return [(ax + t * dx, ay + t * dy)]
    sq = math.sqrt(disc)
    t1 = (-bb + sq) / (2.0 * aa)
    t2 = (-bb - sq) / (2.0 * aa)
    return [(ax + t1 * dx, ay + t1 * dy), (ax + t2 * dx, ay + t2 * dy)]


def circle_circle_intersection(c1, c2):
    """Intersection points of two circles ``(cx, cy, r)``.

    Returns a list of 0, 1 (tangent), or 2 ``(x, y)`` points. Coincident circles raise
    (infinitely many intersections).
    """
    x1, y1, r1 = c1
    x2, y2, r2 = c2
    dx, dy = x2 - x1, y2 - y1
    dcen = math.hypot(dx, dy)
    if dcen == 0.0 and r1 == r2:
        raise ValueError("coincident circles intersect everywhere")
    if dcen > r1 + r2 + 1e-12 or dcen < abs(r1 - r2) - 1e-12:
        return []                                   # separate or one inside the other
    # Distance from c1 to the radical line, and the half-chord height.
    a = (r1 * r1 - r2 * r2 + dcen * dcen) / (2.0 * dcen)
    h2 = r1 * r1 - a * a
    xm = x1 + a * dx / dcen
    ym = y1 + a * dy / dcen
    if h2 <= 1e-12:
        return [(xm, ym)]                           # tangent
    h = math.sqrt(h2)
    rx = -dy * (h / dcen)
    ry = dx * (h / dcen)
    return [(xm + rx, ym + ry), (xm - rx, ym - ry)]
