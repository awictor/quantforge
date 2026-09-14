"""Planar computational geometry: convex hull, area, point-in-polygon, centroid.

The core 2-D geometry primitives on lists of ``(x, y)`` points: the convex hull by
Andrew's monotone chain, signed polygon area and centroid by the shoelace formula,
point-in-polygon by ray casting, and the closest pair of points. These underpin spatial
queries, collision tests, and any "shape of a point cloud" computation. Pure standard
library.
"""

import math


def _cross(o, a, b):
    """Cross product ``(a - o) x (b - o)`` -- positive if o->a->b turns left."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def convex_hull(points):
    """Convex hull of a set of 2-D points (Andrew's monotone chain).

    Returns the hull vertices in counter-clockwise order, starting from the lowest-then-
    leftmost point, without repeating the first point. Collinear interior points are
    dropped. Needs at least one point; duplicates are ignored.
    """
    pts = sorted(set(map(tuple, points)))
    if len(pts) <= 2:
        return pts
    lower = []
    for p in pts:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    # Concatenate, dropping each list's last point (shared with the other).
    return lower[:-1] + upper[:-1]


def polygon_area(polygon):
    """Area of a simple polygon by the shoelace formula (unsigned).

    ``polygon`` is a list of ``(x, y)`` vertices in order (open ring; the last vertex is
    joined back to the first). Returns the absolute area. Needs at least three vertices.
    """
    n = len(polygon)
    if n < 3:
        raise ValueError("polygon needs at least three vertices")
    s = 0.0
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) / 2.0


def polygon_centroid(polygon):
    """Centroid ``(cx, cy)`` of a simple polygon (area-weighted, shoelace form).

    Falls back to the vertex average for a degenerate (zero-area) polygon. ``polygon`` is
    an ordered ``(x, y)`` vertex list.
    """
    n = len(polygon)
    if n < 3:
        raise ValueError("polygon needs at least three vertices")
    a = 0.0
    cx = 0.0
    cy = 0.0
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        cross = x1 * y2 - x2 * y1
        a += cross
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
    if a == 0.0:
        return (sum(p[0] for p in polygon) / n, sum(p[1] for p in polygon) / n)
    a *= 0.5
    return (cx / (6.0 * a), cy / (6.0 * a))


def point_in_polygon(point, polygon):
    """Whether ``point`` lies inside a simple ``polygon`` (ray-casting, odd-crossing rule).

    Casts a ray to the right and counts edge crossings; an odd count means inside. Points
    exactly on an edge are reported as inside. ``polygon`` is an ordered ``(x, y)`` vertex
    list.
    """
    n = len(polygon)
    if n < 3:
        raise ValueError("polygon needs at least three vertices")
    x, y = point
    inside = False
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        # On-edge check (collinear and within the segment's bounding box).
        if _on_segment(point, polygon[i], polygon[(i + 1) % n]):
            return True
        if (y1 > y) != (y2 > y):
            x_cross = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x < x_cross:
                inside = not inside
    return inside


def _on_segment(p, a, b):
    if abs((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])) > 1e-12:
        return False
    return (min(a[0], b[0]) - 1e-12 <= p[0] <= max(a[0], b[0]) + 1e-12 and
            min(a[1], b[1]) - 1e-12 <= p[1] <= max(a[1], b[1]) + 1e-12)


def closest_pair(points):
    """Closest pair of points and their distance: ``(p, q, distance)``.

    Divide-and-conquer in ``O(n log n)``. Needs at least two points; ties break to the
    first pair found.
    """
    n = len(points)
    if n < 2:
        raise ValueError("need at least two points")
    pts = sorted(map(tuple, points))

    def dist(a, b):
        return math.hypot(a[0] - b[0], a[1] - b[1])

    def rec(px):
        m = len(px)
        if m <= 3:
            best = None
            for i in range(m):
                for j in range(i + 1, m):
                    d = dist(px[i], px[j])
                    if best is None or d < best[2]:
                        best = (px[i], px[j], d)
            return best
        mid = m // 2
        midx = px[mid][0]
        left = rec(px[:mid])
        right = rec(px[mid:])
        best = left if left[2] <= right[2] else right
        strip = [p for p in px if abs(p[0] - midx) < best[2]]
        strip.sort(key=lambda p: p[1])
        for i in range(len(strip)):
            for j in range(i + 1, min(i + 7, len(strip))):
                d = dist(strip[i], strip[j])
                if d < best[2]:
                    best = (strip[i], strip[j], d)
        return best

    return rec(pts)
