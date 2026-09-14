"""Delaunay triangulation of a 2-D point set (Bowyer-Watson algorithm).

The Delaunay triangulation is the one in which no point lies inside any triangle's
circumcircle -- it maximizes the minimum angle, avoiding thin slivers, and is the dual of
the Voronoi diagram. Bowyer-Watson builds it incrementally: each new point deletes the
triangles whose circumcircle contains it, leaving a polygonal hole that is re-triangulated
to the point. Returns triangles as index triples into the input list. Pure standard library.
"""


def _circumcircle_contains(a, b, c, p):
    # True if p is strictly inside the circumcircle of triangle (a, b, c).
    # Uses the standard in-circle determinant; assumes (a, b, c) counter-clockwise.
    ax, ay = a[0] - p[0], a[1] - p[1]
    bx, by = b[0] - p[0], b[1] - p[1]
    cx, cy = c[0] - p[0], c[1] - p[1]
    det = (
        (ax * ax + ay * ay) * (bx * cy - cx * by)
        - (bx * bx + by * by) * (ax * cy - cx * ay)
        + (cx * cx + cy * cy) * (ax * by - bx * ay)
    )
    return det > 1e-12


def _ccw(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]) > 0


def delaunay_triangulation(points):
    """Return the Delaunay triangulation of ``points`` as a list of index triples.

    Each triple ``(i, j, k)`` indexes into ``points``. Needs at least three
    non-collinear points. Duplicate points are ignored for the triangulation but indices
    refer to the original list.
    """
    pts = [tuple(p) for p in points]
    n = len(pts)
    if n < 3:
        raise ValueError("need at least three points")

    # super-triangle enclosing all points
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    dx = maxx - minx or 1.0
    dy = maxy - miny or 1.0
    dmax = max(dx, dy)
    midx = (minx + maxx) / 2
    midy = (miny + maxy) / 2
    # three super-vertices, appended after the real points (indices n, n+1, n+2)
    sv = [
        (midx - 20 * dmax, midy - dmax),
        (midx, midy + 20 * dmax),
        (midx + 20 * dmax, midy - dmax),
    ]
    allp = pts + sv
    # each triangle stored as an index triple in CCW order
    def tri(i, j, k):
        a, b, c = allp[i], allp[j], allp[k]
        if not _ccw(a, b, c):
            j, k = k, j
        return (i, j, k)

    tris = [tri(n, n + 1, n + 2)]
    for pi in range(n):
        p = allp[pi]
        bad = [t for t in tris if _circumcircle_contains(allp[t[0]], allp[t[1]], allp[t[2]], p)]
        # boundary of the polygonal hole: edges belonging to exactly one bad triangle
        edge_count = {}
        for (i, j, k) in bad:
            for e in ((i, j), (j, k), (k, i)):
                key = frozenset(e)
                edge_count[key] = edge_count.get(key, (0, None))
                edge_count[key] = (edge_count[key][0] + 1, e)
        boundary = [e for key, (cnt, e) in edge_count.items() if cnt == 1]
        tris = [t for t in tris if t not in bad]
        for (a, b) in boundary:
            tris.append(tri(a, b, pi))
    # drop triangles touching a super-vertex
    result = [t for t in tris if all(idx < n for idx in t)]
    return result
