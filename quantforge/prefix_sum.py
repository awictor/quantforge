"""Prefix sums: O(1) range-sum queries and O(1) range updates on static/offline data.

A prefix-sum table precomputes cumulative sums once, then answers any range-sum query in
constant time by subtracting two entries:

  * ``PrefixSum1D`` -- ``range_sum(lo, hi)`` over a fixed array.
  * ``PrefixSum2D`` -- the summed-area table (integral image): the sum over any axis-aligned
    rectangle by inclusion-exclusion of four corners.
  * ``DifferenceArray`` -- the dual, for the offline "apply many range-add updates, then
    read the final array" pattern: each range add is ``O(1)``, and one final pass
    materializes the result.

Pure standard library.
"""


class PrefixSum1D:
    """Constant-time range sums over a fixed 1-D array.

    ``range_sum(lo, hi)`` returns the sum of ``values[lo:hi]`` (half-open) in ``O(1)`` after
    an ``O(n)`` build. ``prefix(i)`` gives the sum of the first ``i`` elements.
    """

    __slots__ = ("_pre", "_n")

    def __init__(self, values):
        data = list(values)
        self._n = len(data)
        pre = [0] * (self._n + 1)
        for i, v in enumerate(data):
            pre[i + 1] = pre[i] + v
        self._pre = pre

    def __len__(self):
        return self._n

    def prefix(self, i):
        """Sum of the first ``i`` elements (``0 <= i <= n``)."""
        if not 0 <= i <= self._n:
            raise IndexError("prefix index out of range")
        return self._pre[i]

    def range_sum(self, lo, hi):
        """Sum of ``values[lo:hi]`` (half-open, ``0 <= lo <= hi <= n``)."""
        if not 0 <= lo <= hi <= self._n:
            raise IndexError("range out of bounds")
        return self._pre[hi] - self._pre[lo]

    def total(self):
        """Sum of every element."""
        return self._pre[self._n]


class PrefixSum2D:
    """Summed-area table for constant-time rectangle sums over a fixed 2-D grid.

    Build from a list of equal-length rows. ``range_sum(r0, c0, r1, c1)`` returns the sum
    over rows ``[r0, r1)`` and columns ``[c0, c1)`` (half-open) in ``O(1)``.
    """

    __slots__ = ("_pre", "_rows", "_cols")

    def __init__(self, grid):
        rows = len(grid)
        cols = len(grid[0]) if rows else 0
        if any(len(row) != cols for row in grid):
            raise ValueError("all rows must have equal length")
        self._rows = rows
        self._cols = cols
        pre = [[0] * (cols + 1) for _ in range(rows + 1)]
        for i in range(rows):
            row = grid[i]
            for j in range(cols):
                pre[i + 1][j + 1] = (
                    row[j] + pre[i][j + 1] + pre[i + 1][j] - pre[i][j]
                )
        self._pre = pre

    @property
    def shape(self):
        """The ``(rows, cols)`` of the underlying grid."""
        return (self._rows, self._cols)

    def range_sum(self, r0, c0, r1, c1):
        """Sum over rows ``[r0, r1)`` and columns ``[c0, c1)`` (half-open)."""
        if not (0 <= r0 <= r1 <= self._rows and 0 <= c0 <= c1 <= self._cols):
            raise IndexError("rectangle out of bounds")
        p = self._pre
        return p[r1][c1] - p[r0][c1] - p[r1][c0] + p[r0][c0]

    def total(self):
        """Sum of every cell."""
        return self._pre[self._rows][self._cols]


class DifferenceArray:
    """Offline range-add / final-read via a difference array.

    ``add(lo, hi, delta)`` adds ``delta`` to every index in ``[lo, hi)`` in ``O(1)``; after
    all updates, ``result()`` materializes the final array in ``O(n)``. Ideal when many
    range updates precede a single read.
    """

    __slots__ = ("_diff", "_n")

    def __init__(self, n_or_values):
        if isinstance(n_or_values, int):
            self._n = n_or_values
            if self._n < 0:
                raise ValueError("length must be non-negative")
            self._diff = [0] * (self._n + 1)
        else:
            data = list(n_or_values)
            self._n = len(data)
            self._diff = [0] * (self._n + 1)
            for i, v in enumerate(data):
                self._diff[i] += v
                self._diff[i + 1] -= v

    def __len__(self):
        return self._n

    def add(self, lo, hi, delta):
        """Add ``delta`` to every index in ``[lo, hi)`` (half-open) in ``O(1)``."""
        if not 0 <= lo <= hi <= self._n:
            raise IndexError("range out of bounds")
        self._diff[lo] += delta
        self._diff[hi] -= delta
        return self

    def result(self):
        """Materialize the final array after all range-adds (``O(n)``)."""
        out = [0] * self._n
        running = 0
        for i in range(self._n):
            running += self._diff[i]
            out[i] = running
        return out
