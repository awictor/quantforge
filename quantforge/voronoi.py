"""Voronoi diagram from the Delaunay dual: vertices, adjacency, nearest-site.

The Voronoi diagram partitions the plane into cells, one per site, where each cell is the
region closer to its site than to any other. It is the geometric dual of the Delaunay
triangulation: every Delaunay triangle's circumcenter is a Voronoi vertex, and two sites
share a Voronoi edge exactly when they share a Delaunay edge. This module derives the
Voronoi vertices and the site-adjacency (Delaunay neighbour) graph from
:func:`quantforge.delaunay.delaunay_triangulation`, and answers nearest-site queries
directly. Pure standard library.
"""

from .delaunay import delaunay_triangulation


def _circumcenter(a, b, c):
    ax, ay = a
    bx, by = b
    cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-15:
        raise ValueError("collinear triangle has no circumcenter")
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay)
          + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx)
          + (cx * cx + cy * cy) * (bx - ax)) / d
    return (ux, uy)


def voronoi_vertices(points):
    """Return the Voronoi vertices: the circumcenter of each Delaunay triangle.

    Each vertex is equidistant from the three sites of its triangle. Returns a list of
    ``(x, y)`` points (one per Delaunay triangle).
    """
    tris = delaunay_triangulation(points)
    pts = [tuple(p) for p in points]
    return [_circumcenter(pts[i], pts[j], pts[k]) for (i, j, k) in tris]


def delaunay_neighbors(points):
    """Return the site-adjacency graph ``{i: set(neighbor indices)}`` (Delaunay edges).

    Two sites are neighbours iff they share a Delaunay edge -- equivalently, iff their
    Voronoi cells share an edge. The relation is symmetric.
    """
    tris = delaunay_triangulation(points)
    adj = {i: set() for i in range(len(points))}
    for (i, j, k) in tris:
        adj[i].add(j)
        adj[j].add(i)
        adj[j].add(k)
        adj[k].add(j)
        adj[i].add(k)
        adj[k].add(i)
    return adj


def nearest_site(points, query):
    """Return the index of the site nearest to ``query`` (the Voronoi cell it falls in)."""
    if not points:
        raise ValueError("need at least one site")
    qx, qy = query
    best = 0
    best_d = float("inf")
    for idx, (px, py) in enumerate(points):
        d = (px - qx) ** 2 + (py - qy) ** 2
        if d < best_d:
            best_d = d
            best = idx
    return best
