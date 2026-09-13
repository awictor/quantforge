"""QR decomposition by Householder reflections, and QR least squares.

Factors an ``m x n`` matrix (``m >= n``) as ``A = Q R`` with ``Q`` orthogonal
(``Q' Q = I``) and ``R`` upper triangular, using numerically stable Householder
reflections. Solving least squares through QR (``R x = Q' b``) avoids forming the
ill-conditioned normal-equations matrix ``A' A``, so it is the preferred solver for
regression when the design is nearly collinear. Pure standard library.
"""


def qr_decomposition(A):
    """Householder QR of an ``m x n`` matrix ``A`` (``m >= n``).

    Returns ``(Q, R)`` with ``Q`` an ``m x m`` orthogonal matrix and ``R`` an
    ``m x n`` upper-triangular matrix such that ``A = Q R``. Pure Python lists.
    """
    m = len(A)
    if m == 0:
        raise ValueError("A must be non-empty")
    n = len(A[0])
    if any(len(row) != n for row in A):
        raise ValueError("A must be rectangular")
    if m < n:
        raise ValueError("need m >= n rows")

    R = [row[:] for row in A]
    Q = [[1.0 if i == j else 0.0 for j in range(m)] for i in range(m)]

    for k in range(min(m - 1, n)):
        # Householder vector for column k below the diagonal.
        x = [R[i][k] for i in range(k, m)]
        norm_x = sum(v * v for v in x) ** 0.5
        if norm_x == 0.0:
            continue
        alpha = -norm_x if x[0] >= 0 else norm_x
        v = x[:]
        v[0] -= alpha
        vnorm = sum(t * t for t in v) ** 0.5
        if vnorm == 0.0:
            continue
        v = [t / vnorm for t in v]

        # Apply H = I - 2 v v' to R (rows k..m-1).
        for j in range(n):
            dot = sum(v[i] * R[k + i][j] for i in range(m - k))
            for i in range(m - k):
                R[k + i][j] -= 2.0 * v[i] * dot
        # Accumulate Q (apply H on the right of Q, columns k..m-1).
        for i in range(m):
            dot = sum(Q[i][k + t] * v[t] for t in range(m - k))
            for t in range(m - k):
                Q[i][k + t] -= 2.0 * dot * v[t]

    # Clean tiny sub-diagonal noise in R.
    for i in range(m):
        for j in range(min(i, n)):
            R[i][j] = 0.0
    return Q, R


def qr_solve(A, b):
    """Least-squares solution of ``A x = b`` via QR (``x`` minimizes ||A x - b||).

    Computes ``Q' b`` and back-substitutes the upper-triangular ``R``. Numerically
    stabler than the normal equations for ill-conditioned ``A``. ``A`` is ``m x n``
    with ``m >= n``; returns the length-``n`` coefficient vector.
    """
    m = len(A)
    n = len(A[0])
    if len(b) != m:
        raise ValueError("b length must match the rows of A")
    Q, R = qr_decomposition(A)
    # c = Q' b (first n entries suffice).
    c = [sum(Q[i][j] * b[i] for i in range(m)) for j in range(n)]
    # Back-substitution on the n x n upper-triangular top block of R.
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        if R[i][i] == 0.0:
            raise ValueError("A is rank-deficient; R has a zero pivot")
        s = c[i] - sum(R[i][j] * x[j] for j in range(i + 1, n))
        x[i] = s / R[i][i]
    return x
