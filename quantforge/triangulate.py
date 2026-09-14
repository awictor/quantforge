"""Polygon triangulation by ear clipping, plus orientation and convexity tests.

Any simple polygon can be cut into ``n - 2`` triangles; ear clipping does it in
``O(n^2)`` by repeatedly snipping off "ears" (a vertex whose triangle lies inside the
polygon and contains no other vertex). Triangulation underpins area computation,
rendering, and finite-element meshing. Also provides orientation (clockwise vs
counter-clockwise) and a convexity test. Points are ``(x, y)`` tuples. Pure standard
library.
"""


def signed_area(polygon):
    """Signed area of a polygon (shoelace): positive if counter-clockwise, negative if CW."""
    n = len(polygon)
    if n < 3:
        raise ValueError("polygon needs at least three vertices")
    s = 0.0
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return s / 2.0


def is_clockwise(polygon):
    """Whether the polygon vertices are ordered clockwise (negative signed area)."""
    return signed_area(polygon) < 0.0


def is_convex_polygon(polygon):
    """Whether a simple polygon is convex (all turns the same direction)."""
    n = len(polygon)
    if n < 3:
        raise ValueError("polygon needs at least three vertices")
    sign = 0
    for i in range(n):
        ax, ay = polygon[i]
        bx, by = polygon[(i + 1) % n]
        cx, cy = polygon[(i + 2) % n]
        cross = (bx - ax) * (cy - by) - (by - ay) * (cx - bx)
        if cross != 0:
            s = 1 if cross > 0 else -1
            if sign == 0:
                sign = s
            elif s != sign:
                return False
    return True


def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _point_in_triangle(p, a, b, c):
    d1 = _cross(a, b, p)
    d2 = _cross(b, c, p)
    d3 = _cross(c, a, p)
    has_neg = d1 < 0 or d2 < 0 or d3 < 0
    has_pos = d1 > 0 or d2 > 0 or d3 > 0
    return not (has_neg and has_pos)


def ear_clipping_triangulate(polygon):
    """Triangulate a simple polygon by ear clipping.

    Returns a list of ``n - 2`` triangles, each a tuple of three ``(x, y)`` vertices. The
    polygon must be simple (non-self-intersecting) with vertices in order; both windings
    are accepted. The triangles' areas sum to the polygon's area.
    """
    n = len(polygon)
    if n < 3:
        raise ValueError("polygon needs at least three vertices")
    # Work with a CCW copy so the interior is consistently to the left.
    pts = list(polygon)
    if is_clockwise(pts):
        pts = pts[::-1]
    idx = list(range(len(pts)))
    triangles = []
    guard = 0
    max_guard = len(idx) ** 2 + 10
    while len(idx) > 3:
        guard += 1
        if guard > max_guard:
            raise ValueError("triangulation failed (polygon may be non-simple)")
        ear_found = False
        m = len(idx)
        for i in range(m):
            ia, ib, ic = idx[(i - 1) % m], idx[i], idx[(i + 1) % m]
            a, b, c = pts[ia], pts[ib], pts[ic]
            if _cross(a, b, c) <= 0:
                continue                       # reflex or collinear vertex, not an ear tip
            # No other polygon vertex inside triangle a-b-c.
            if any(_point_in_triangle(pts[j], a, b, c)
                   for j in idx if j not in (ia, ib, ic)):
                continue
            triangles.append((a, b, c))
            del idx[i]
            ear_found = True
            break
        if not ear_found:
            raise ValueError("triangulation failed (polygon may be non-simple)")
    triangles.append((pts[idx[0]], pts[idx[1]], pts[idx[2]]))
    return triangles
