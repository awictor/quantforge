"""2-D Fenwick tree (binary indexed tree): dynamic point updates and rectangle sums.

Where `PrefixSum2D` gives constant-time rectangle sums on a *static* grid, the 2-D Fenwick
tree supports *updates*: add to a cell and query the sum over any axis-aligned rectangle,
each in ``O(log rows * log cols)``. It generalizes the 1-D binary indexed tree to two
dimensions by nesting the ``i & -i`` traversal. Ideal for running 2-D counts/sums that
change over time (collision grids, dynamic heatmaps). Pure standard library.
"""


class FenwickTree2D:
    """Dynamic 2-D prefix sums: point-add and rectangle-sum in ``O(log R * log C)``.

    Construct with the grid shape ``(rows, cols)`` (optionally an initial grid). ``add(r, c,
    delta)`` adds to a cell; ``prefix_sum(r, c)`` sums the rectangle ``[0, r) x [0, c)``; and
    ``range_sum(r0, c0, r1, c1)`` sums ``[r0, r1) x [c0, c1)`` by inclusion-exclusion.
    """

    __slots__ = ("_rows", "_cols", "_tree")

    def __init__(self, rows=None, cols=None, grid=None):
        if grid is not None:
            rows = len(grid)
            cols = len(grid[0]) if rows else 0
            if any(len(row) != cols for row in grid):
                raise ValueError("all rows must have equal length")
        if rows is None or cols is None:
            raise ValueError("rows and cols required when no grid is given")
        if rows < 0 or cols < 0:
            raise ValueError("dimensions must be non-negative")
        self._rows = rows
        self._cols = cols
        self._tree = [[0] * (cols + 1) for _ in range(rows + 1)]
        if grid is not None:
            for i in range(rows):
                for j in range(cols):
                    if grid[i][j]:
                        self.add(i, j, grid[i][j])

    @property
    def shape(self):
        """The ``(rows, cols)`` of the grid."""
        return (self._rows, self._cols)

    def add(self, r, c, delta):
        """Add ``delta`` to cell ``(r, c)``."""
        if not (0 <= r < self._rows and 0 <= c < self._cols):
            raise IndexError("cell out of bounds")
        i = r + 1
        while i <= self._rows:
            j = c + 1
            while j <= self._cols:
                self._tree[i][j] += delta
                j += j & -j
            i += i & -i

    def prefix_sum(self, r, c):
        """Sum over the rectangle ``[0, r) x [0, c)`` (half-open)."""
        if not (0 <= r <= self._rows and 0 <= c <= self._cols):
            raise IndexError("prefix bounds out of range")
        total = 0
        i = r
        while i > 0:
            j = c
            while j > 0:
                total += self._tree[i][j]
                j -= j & -j
            i -= i & -i
        return total

    def range_sum(self, r0, c0, r1, c1):
        """Sum over ``[r0, r1) x [c0, c1)`` (half-open) by four-corner inclusion-exclusion."""
        if not (0 <= r0 <= r1 <= self._rows and 0 <= c0 <= c1 <= self._cols):
            raise IndexError("rectangle out of bounds")
        return (
            self.prefix_sum(r1, c1)
            - self.prefix_sum(r0, c1)
            - self.prefix_sum(r1, c0)
            + self.prefix_sum(r0, c0)
        )
