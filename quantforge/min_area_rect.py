"""Minimum-area enclosing rectangle of a 2-D point set (rotating calipers).

The smallest-area rectangle containing a set of points has one side flush with an edge of
their convex hull (a theorem of Freeman & Shapira), so it suffices to test the rectangle
aligned to each hull edge and keep the smallest. Projecting the hull onto each edge's
direction and its normal gives that rectangle's extent in ``O(h)`` per edge, ``O(h^2)``
overall for ``h`` hull vertices. Returns the area, side lengths, orientation, and corners.
Pure standard library.
"""

import math

from .geometry import convex_hull


def min_area_rectangle(points):
    """Return the minimum-area enclosing rectangle of ``points``.

    Result is a dict with ``area``, ``width``, ``height`` (the two side lengths, width the
    longer), ``angle`` (edge direction in radians), and ``corners`` (four ``(x, y)`` points
    in order). Needs at least one point; a degenerate (collinear) set gives a zero-area
    rectangle.
    """
    hull = convex_hull(points)
    if not hull:
        raise ValueError("need at least one point")
    if len(hull) == 1:
        p = hull[0]
        return {"area": 0.0, "width": 0.0, "height": 0.0, "angle": 0.0,
                "corners": [p, p, p, p]}
    if len(hull) == 2:
        (x0, y0), (x1, y1) = hull
        length = math.hypot(x1 - x0, y1 - y0)
        angle = math.atan2(y1 - y0, x1 - x0)
        return {"area": 0.0, "width": length, "height": 0.0, "angle": angle,
                "corners": [hull[0], hull[1], hull[1], hull[0]]}
    n = len(hull)
    best = None
    for i in range(n):
        ax, ay = hull[i]
        bx, by = hull[(i + 1) % n]
        ex, ey = bx - ax, by - ay
        elen = math.hypot(ex, ey)
        if elen == 0:
            continue
        ux, uy = ex / elen, ey / elen        # edge unit direction
        nx, ny = -uy, ux                     # normal
        min_u = min_n = float("inf")
        max_u = max_n = float("-inf")
        for px, py in hull:
            du = (px - ax) * ux + (py - ay) * uy
            dn = (px - ax) * nx + (py - ay) * ny
            min_u, max_u = min(min_u, du), max(max_u, du)
            min_n, max_n = min(min_n, dn), max(max_n, dn)
        w = max_u - min_u
        h = max_n - min_n
        area = w * h
        if best is None or area < best[0]:
            corners = [
                (ax + min_u * ux + min_n * nx, ay + min_u * uy + min_n * ny),
                (ax + max_u * ux + min_n * nx, ay + max_u * uy + min_n * ny),
                (ax + max_u * ux + max_n * nx, ay + max_u * uy + max_n * ny),
                (ax + min_u * ux + max_n * nx, ay + min_u * uy + max_n * ny),
            ]
            best = (area, w, h, math.atan2(uy, ux), corners)
    area, w, h, angle, corners = best
    width, height = (w, h) if w >= h else (h, w)
    return {"area": area, "width": width, "height": height, "angle": angle,
            "corners": corners}
