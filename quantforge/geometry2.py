"""Planar geometry II: segment intersection and convex polygon clipping.

Builds on :mod:`quantforge.geometry` with the operations behind collision detection and
map/window clipping: whether two segments cross (and where), the perimeter of a polygon,
and Sutherland-Hodgman clipping of a polygon against a convex window. Pure standard
library.
"""

import math


def _orient(a, b, c):
    """Orientation of the ordered triple: >0 left turn, <0 right turn, 0 collinear."""
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _on_segment(a, b, p):
    """Whether collinear point ``p`` lies on segment ``a-b``."""
    return (min(a[0], b[0]) - 1e-12 <= p[0] <= max(a[0], b[0]) + 1e-12 and
            min(a[1], b[1]) - 1e-12 <= p[1] <= max(a[1], b[1]) + 1e-12)


def segments_intersect(p1, p2, p3, p4):
    """Whether segment ``p1-p2`` intersects segment ``p3-p4`` (including endpoints/collinear).

    Uses the standard four-orientation test with collinear-overlap handling. Returns a
    bool.
    """
    d1 = _orient(p3, p4, p1)
    d2 = _orient(p3, p4, p2)
    d3 = _orient(p1, p2, p3)
    d4 = _orient(p1, p2, p4)
    if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)):
        # proper crossing (strict opposite sides both ways)
        if d1 != 0 and d2 != 0 and d3 != 0 and d4 != 0:
            return True
    # Collinear / touching cases.
    if d1 == 0 and _on_segment(p3, p4, p1):
        return True
    if d2 == 0 and _on_segment(p3, p4, p2):
        return True
    if d3 == 0 and _on_segment(p1, p2, p3):
        return True
    if d4 == 0 and _on_segment(p1, p2, p4):
        return True
    return False


def segment_intersection(p1, p2, p3, p4):
    """Intersection point of two segments, or ``None`` if they do not cross at a single point.

    Returns the ``(x, y)`` crossing for segments that meet at exactly one point; returns
    ``None`` if they are parallel, collinear, or disjoint. Endpoints count as
    intersections.
    """
    r = (p2[0] - p1[0], p2[1] - p1[1])
    s = (p4[0] - p3[0], p4[1] - p3[1])
    denom = r[0] * s[1] - r[1] * s[0]
    if denom == 0.0:
        return None                      # parallel or collinear
    qp = (p3[0] - p1[0], p3[1] - p1[1])
    t = (qp[0] * s[1] - qp[1] * s[0]) / denom
    u = (qp[0] * r[1] - qp[1] * r[0]) / denom
    if -1e-12 <= t <= 1 + 1e-12 and -1e-12 <= u <= 1 + 1e-12:
        return (p1[0] + t * r[0], p1[1] + t * r[1])
    return None


def polygon_perimeter(polygon):
    """Perimeter of a polygon: the sum of its edge lengths (closed ring)."""
    n = len(polygon)
    if n < 2:
        raise ValueError("polygon needs at least two vertices")
    total = 0.0
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        total += math.hypot(x2 - x1, y2 - y1)
    return total


def clip_polygon(subject, clip):
    """Clip ``subject`` polygon against a convex ``clip`` polygon (Sutherland-Hodgman).

    Both are ordered ``(x, y)`` vertex lists; ``clip`` must be convex and given
    counter-clockwise. Returns the clipped polygon (possibly empty) as a vertex list --
    the intersection of the subject with the clip window. The classic viewport/window
    clipping algorithm.
    """
    if len(clip) < 3:
        raise ValueError("clip polygon needs at least three vertices")
    output = list(subject)
    cn = len(clip)
    for i in range(cn):
        a = clip[i]
        b = clip[(i + 1) % cn]
        if not output:
            break
        inp = output
        output = []
        for j in range(len(inp)):
            cur = inp[j]
            prev = inp[j - 1]
            cur_in = _orient(a, b, cur) >= 0
            prev_in = _orient(a, b, prev) >= 0
            if cur_in:
                if not prev_in:
                    output.append(_line_intersection(prev, cur, a, b))
                output.append(cur)
            elif prev_in:
                output.append(_line_intersection(prev, cur, a, b))
    return output


def _line_intersection(p1, p2, a, b):
    """Intersection of segment ``p1-p2`` with the *infinite line* through ``a`` and ``b``.

    Used by Sutherland-Hodgman, which clips each subject edge against the full line of a
    clip edge (not just the clip segment). ``p1-p2`` is known to straddle the line.
    """
    r = (p2[0] - p1[0], p2[1] - p1[1])
    s = (b[0] - a[0], b[1] - a[1])
    denom = r[0] * s[1] - r[1] * s[0]
    if denom == 0.0:
        return p2                        # parallel (shouldn't happen for a straddling edge)
    qp = (a[0] - p1[0], a[1] - p1[1])
    t = (qp[0] * s[1] - qp[1] * s[0]) / denom
    return (p1[0] + t * r[0], p1[1] + t * r[1])
