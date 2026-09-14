"""KD-tree: nearest-neighbour and range queries over k-dimensional points.

A k-d tree recursively partitions points by alternating coordinate axes, so a
nearest-neighbour or radius query prunes whole subtrees whose bounding slab cannot hold a
closer point -- typically ``O(log n)`` per query on well-distributed data, versus the
``O(n)`` scan of a brute-force search. Supports ``k``-nearest, all points within a radius,
and axis-aligned box queries. Distances are Euclidean. Pure standard library.
"""

import heapq


def _dist_sq(a, b):
    return sum((x - y) * (x - y) for x, y in zip(a, b))


class _Node:
    __slots__ = ("point", "index", "axis", "left", "right")

    def __init__(self, point, index, axis):
        self.point = point
        self.index = index
        self.axis = axis
        self.left = None
        self.right = None


class KDTree:
    """Static k-d tree over a fixed set of points for spatial queries.

    Build once from a list of equal-length coordinate tuples; then run `nearest`,
    `k_nearest`, `within_radius`, and `range_search`. Query methods return the original
    point *indices* (in build order), so callers can map back to their own payloads.
    """

    __slots__ = ("_root", "_dim", "_n", "_points")

    def __init__(self, points):
        self._points = [tuple(float(c) for c in p) for p in points]
        self._n = len(self._points)
        if self._n == 0:
            self._root = None
            self._dim = 0
            return
        self._dim = len(self._points[0])
        if any(len(p) != self._dim for p in self._points):
            raise ValueError("all points must have the same dimension")
        self._root = self._build(list(range(self._n)), 0)

    def _build(self, idxs, depth):
        if not idxs:
            return None
        axis = depth % self._dim
        idxs.sort(key=lambda i: self._points[i][axis])
        mid = len(idxs) // 2
        node = _Node(self._points[idxs[mid]], idxs[mid], axis)
        node.left = self._build(idxs[:mid], depth + 1)
        node.right = self._build(idxs[mid + 1:], depth + 1)
        return node

    def __len__(self):
        return self._n

    def nearest(self, target):
        """Return ``(index, distance)`` of the single closest point to ``target``."""
        res = self.k_nearest(target, 1)
        if not res:
            raise ValueError("tree is empty")
        return res[0]

    def k_nearest(self, target, k):
        """Return the ``k`` closest points as ``(index, distance)``, ascending by distance."""
        if k < 0:
            raise ValueError("k must be non-negative")
        target = tuple(float(c) for c in target)
        if self._root is None or k == 0:
            return []
        # max-heap of (-dist_sq, index) holding the current k best
        heap = []

        def visit(node):
            if node is None:
                return
            d = _dist_sq(target, node.point)
            if len(heap) < k:
                heapq.heappush(heap, (-d, node.index))
            elif d < -heap[0][0]:
                heapq.heapreplace(heap, (-d, node.index))
            axis = node.axis
            diff = target[axis] - node.point[axis]
            near, far = (node.left, node.right) if diff < 0 else (node.right, node.left)
            visit(near)
            # only descend the far side if its slab could hold a closer point
            if len(heap) < k or diff * diff < -heap[0][0]:
                visit(far)

        visit(self._root)
        out = [(idx, d_sq ** 0.5) for d_sq, idx in ((-nd, i) for nd, i in heap)]
        out.sort(key=lambda t: (t[1], t[0]))
        return out

    def within_radius(self, target, radius):
        """Return ``(index, distance)`` for all points within ``radius`` (inclusive)."""
        if radius < 0:
            raise ValueError("radius must be non-negative")
        target = tuple(float(c) for c in target)
        r_sq = radius * radius
        found = []

        def visit(node):
            if node is None:
                return
            d = _dist_sq(target, node.point)
            if d <= r_sq:
                found.append((node.index, d ** 0.5))
            axis = node.axis
            diff = target[axis] - node.point[axis]
            near, far = (node.left, node.right) if diff < 0 else (node.right, node.left)
            visit(near)
            if diff * diff <= r_sq:
                visit(far)

        visit(self._root)
        found.sort(key=lambda t: (t[1], t[0]))
        return found

    def range_search(self, lower, upper):
        """Return indices of points inside the axis-aligned box ``[lower, upper]``.

        ``lower`` and ``upper`` are per-axis bounds (inclusive). Results are sorted by index.
        """
        if len(lower) != self._dim or len(upper) != self._dim:
            raise ValueError("bounds must match the tree dimension")
        lower = tuple(float(c) for c in lower)
        upper = tuple(float(c) for c in upper)
        found = []

        def visit(node):
            if node is None:
                return
            p = node.point
            if all(lower[a] <= p[a] <= upper[a] for a in range(self._dim)):
                found.append(node.index)
            axis = node.axis
            if lower[axis] <= p[axis]:      # left subtree may contain in-range points
                visit(node.left)
            if p[axis] <= upper[axis]:      # right subtree may contain in-range points
                visit(node.right)

        visit(self._root)
        found.sort()
        return found
