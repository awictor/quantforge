"""Matrix exponential by scaling-and-squaring with a Pade approximant.

``exp(A) = sum_k A^k / k!`` is computed stably by the standard scaling-and-squaring
method: scale ``A`` by ``2^{-s}`` so its norm is small, evaluate a Pade(6,6)
rational approximant of the exponential on the scaled matrix, then square the result
``s`` times. The main use here is turning a continuous-time Markov generator ``Q``
into the transition matrix ``P(t) = exp(Q t)`` (e.g. credit-rating migration).
Pure standard library.
"""


def _matmul(A, B):
    n = len(A)
    m = len(B[0])
    k = len(B)
    return [[sum(A[i][p] * B[p][j] for p in range(k)) for j in range(m)]
            for i in range(n)]


def _identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def _add(A, B, sb=1.0):
    n = len(A)
    return [[A[i][j] + sb * B[i][j] for j in range(len(A[0]))] for i in range(n)]


def _scale(A, c):
    return [[c * v for v in row] for row in A]


def _inf_norm(A):
    return max(sum(abs(v) for v in row) for row in A)


def _solve_matrix(A, B):
    """Solve A X = B for X via Gauss-Jordan (A square, B same row count)."""
    from .portopt import _invert
    return _matmul(_invert(A), B)


def matrix_exp(A):
    """Matrix exponential ``exp(A)`` by scaling-and-squaring with Pade(6,6).

    ``A`` is a square matrix. Returns ``exp(A)``; ``exp(0) = I``, ``exp`` of a
    diagonal matrix is the diagonal of exponentials, and it satisfies the defining
    series. Accurate across a wide norm range thanks to the scaling step.
    """
    n = len(A)
    if n == 0 or any(len(row) != n for row in A):
        raise ValueError("A must be a non-empty square matrix")

    # Scaling: bring the norm below 0.5.
    norm = _inf_norm(A)
    s = 0
    while norm > 0.5:
        norm *= 0.5
        s += 1
    M = _scale(A, 0.5 ** s)

    # Pade(6,6) coefficients c_k = (2p-k)! p! / ((2p)! k! (p-k)!) for p = 6.
    c = [1.0, 1.0 / 2, 5.0 / 44, 1.0 / 66, 1.0 / 792, 1.0 / 15840, 1.0 / 665280]

    # Build U = sum c_odd M^k, V = sum c_even M^k via powers.
    powers = [_identity(n), M]
    for k in range(2, 7):
        powers.append(_matmul(powers[-1], M))
    U = _scale(_identity(n), 0.0)
    V = _scale(_identity(n), 0.0)
    for k in range(7):
        term = _scale(powers[k], c[k])
        if k % 2 == 0:
            V = _add(V, term)
        else:
            U = _add(U, term)
    # exp(M) ~ (V - U)^{-1} (V + U).
    num = _add(V, U)            # V + U
    den = _add(V, U, sb=-1.0)   # V - U
    R = _solve_matrix(den, num)

    # Squaring: R -> R^{2^s}.
    for _ in range(s):
        R = _matmul(R, R)
    return R
