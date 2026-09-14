"""Li Chao tree: the lower (or upper) envelope of a set of lines.

Given many lines ``y = m x + b``, the Li Chao tree answers "what is the minimum ``y`` at this
``x``?" over a fixed set of query abscissae, and lets lines be added in any order -- unlike
the classic convex-hull trick, which needs sorted slopes. Each insert and query is
``O(log n)`` in the number of candidate ``x`` positions. This is the workhorse for the
DP-optimization pattern ``dp[i] = min_j (m_j x_i + b_j)``. Set ``maximize=True`` for the
upper envelope. Pure standard library.
"""


class LiChaoTree:
    """Lower/upper envelope of lines over a fixed sorted list of query points.

    Construct with the ``xs`` at which queries will be made (any reals, deduplicated and
    sorted internally). ``add_line(m, b)`` inserts ``y = m x + b``; ``query(x)`` returns the
    minimum (or maximum, if ``maximize``) ``y`` over all inserted lines at that ``x``.
    ``x`` must be one of the construction points.
    """

    __slots__ = ("_xs", "_pos", "_n", "_lines", "_maximize")

    _INF_LINE = None  # sentinel for an empty node

    def __init__(self, xs, maximize=False):
        pts = sorted(set(xs))
        if not pts:
            raise ValueError("need at least one query point")
        self._xs = pts
        self._pos = {x: i for i, x in enumerate(pts)}
        self._n = len(pts)
        self._lines = [None] * (4 * self._n)
        self._maximize = maximize

    def _better(self, a, b):
        # is value a better than b under the current objective?
        return a > b if self._maximize else a < b

    def _eval(self, line, x):
        return line[0] * x + line[1]

    def add_line(self, m, b):
        """Insert the line ``y = m x + b`` into the envelope."""
        self._insert(1, 0, self._n - 1, (m, b))

    def _insert(self, node, lo, hi, line):
        cur = self._lines[node]
        if cur is None:
            self._lines[node] = line
            return
        mid = (lo + hi) // 2
        x_lo, x_mid = self._xs[lo], self._xs[mid]
        # decide which line wins at the midpoint; keep it, recurse with the loser
        if self._better(self._eval(line, x_mid), self._eval(cur, x_mid)):
            self._lines[node], line, cur = line, cur, line
        if lo == hi:
            return
        # if the new (now-loser) line wins at an endpoint, it dominates that half
        if self._better(self._eval(line, x_lo), self._eval(self._lines[node], x_lo)):
            self._insert(2 * node, lo, mid, line)
        else:
            self._insert(2 * node + 1, mid + 1, hi, line)

    def query(self, x):
        """Return the best (min, or max if ``maximize``) ``y`` over all lines at ``x``.

        ``x`` must be one of the points supplied at construction.
        """
        if x not in self._pos:
            raise ValueError("x was not among the construction points")
        idx = self._pos[x]
        return self._query(1, 0, self._n - 1, idx, x)

    def _query(self, node, lo, hi, idx, x):
        line = self._lines[node]
        best = None if line is None else self._eval(line, x)
        if lo == hi:
            return best
        mid = (lo + hi) // 2
        if idx <= mid:
            child = self._query(2 * node, lo, mid, idx, x)
        else:
            child = self._query(2 * node + 1, mid + 1, hi, idx, x)
        if child is None:
            return best
        if best is None:
            return child
        return child if self._better(child, best) else best
