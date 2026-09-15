"""Compressed sparse row (CSR) matrix storage and operations.

Most large matrices in scientific computing -- finite-element stiffness matrices, graph
adjacency, discretized operators -- are *sparse*: almost all entries are zero. Storing them
densely wastes ``O(n^2)`` memory and makes a matrix-vector product ``O(n^2)`` when it should be
``O(nnz)``. CSR (compressed sparse row) stores only the nonzeros as three arrays -- ``data``, the
column ``indices``, and the row ``indptr`` -- giving ``O(nnz)`` storage and matvec, which is
exactly the operation the Krylov solvers (:mod:`quantforge.gmres`, :mod:`quantforge.bicgstab`,
:mod:`quantforge.lanczos`) need. This makes those matrix-free solvers usable on genuinely large
sparse systems.

:class:`CSRMatrix` builds from a dense matrix or ``(row, col, value)`` triplets and supports
matrix-vector product, transpose, and conversion back to dense. Pure standard library.
"""


class CSRMatrix:
    """A matrix in compressed-sparse-row form.

    Construct with :meth:`from_dense` or :meth:`from_triplets`. ``matvec`` is ``O(nnz)`` and is
    directly usable as the operator argument to the Krylov solvers.
    """

    def __init__(self, data, indices, indptr, shape):
        self.data = data          # nonzero values, row by row
        self.indices = indices    # column index of each value
        self.indptr = indptr      # row r occupies data[indptr[r]:indptr[r+1]]
        self.shape = shape        # (nrows, ncols)

    @classmethod
    def from_dense(cls, A, tol=0.0):
        """Build a CSR matrix from a dense list-of-rows, dropping entries with ``|a| <= tol``."""
        m = len(A)
        n = len(A[0]) if m else 0
        data, indices, indptr = [], [], [0]
        for i in range(m):
            for j in range(n):
                if abs(A[i][j]) > tol:
                    data.append(A[i][j])
                    indices.append(j)
            indptr.append(len(data))
        return cls(data, indices, indptr, (m, n))

    @classmethod
    def from_triplets(cls, triplets, shape):
        """Build from ``(row, col, value)`` triplets; duplicate ``(row, col)`` pairs are summed."""
        m, n = shape
        rows = [{} for _ in range(m)]
        for r, c, v in triplets:
            if not (0 <= r < m and 0 <= c < n):
                raise ValueError("triplet index out of range")
            rows[r][c] = rows[r].get(c, 0.0) + v
        data, indices, indptr = [], [], [0]
        for i in range(m):
            for c in sorted(rows[i]):
                if rows[i][c] != 0.0:
                    data.append(rows[i][c])
                    indices.append(c)
            indptr.append(len(data))
        return cls(data, indices, indptr, shape)

    @property
    def nnz(self):
        """Number of stored (nonzero) entries."""
        return len(self.data)

    def matvec(self, x):
        """Matrix-vector product ``A @ x`` in ``O(nnz)``."""
        m, n = self.shape
        if len(x) != n:
            raise ValueError("dimension mismatch")
        out = [0.0] * m
        for i in range(m):
            s = 0.0
            for k in range(self.indptr[i], self.indptr[i + 1]):
                s += self.data[k] * x[self.indices[k]]
            out[i] = s
        return out

    def __matmul__(self, x):
        return self.matvec(x)

    def transpose(self):
        """Return the transpose as a new :class:`CSRMatrix` (CSC-of-A built as CSR-of-A^T)."""
        m, n = self.shape
        # count entries per column (= rows of the transpose)
        counts = [0] * n
        for c in self.indices:
            counts[c] += 1
        indptr = [0] * (n + 1)
        for c in range(n):
            indptr[c + 1] = indptr[c] + counts[c]
        data = [0.0] * self.nnz
        indices = [0] * self.nnz
        nextpos = indptr[:n]
        for i in range(m):
            for k in range(self.indptr[i], self.indptr[i + 1]):
                c = self.indices[k]
                p = nextpos[c]
                data[p] = self.data[k]
                indices[p] = i
                nextpos[c] += 1
        return CSRMatrix(data, indices, indptr, (n, m))

    def to_dense(self):
        """Reconstruct the dense list-of-rows matrix."""
        m, n = self.shape
        A = [[0.0] * n for _ in range(m)]
        for i in range(m):
            for k in range(self.indptr[i], self.indptr[i + 1]):
                A[i][self.indices[k]] = self.data[k]
        return A
