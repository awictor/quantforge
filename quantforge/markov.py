"""Finite-state Markov chains: transitions, stationary law, hitting times.

Discrete-time Markov chains over a finite state space, given a row-stochastic
transition matrix ``P`` (row ``i`` = distribution of the next state from ``i``).
Provides the n-step transition matrix, the stationary distribution (the left
eigenvector of ``P`` for eigenvalue 1, found by power iteration), and expected
first-passage / hitting times. Pure standard library.
"""


def _validate_matrix(P):
    n = len(P)
    if n == 0:
        raise ValueError("transition matrix must be non-empty")
    for row in P:
        if len(row) != n:
            raise ValueError("transition matrix must be square")
        if any(x < 0 for x in row):
            raise ValueError("transition probabilities must be non-negative")
        if abs(sum(row) - 1.0) > 1e-9:
            raise ValueError("each row must sum to 1")
    return n


def _matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)]
            for i in range(n)]


def n_step_transition(P, n):
    """``n``-step transition matrix ``P^n`` by repeated squaring.

    Row ``i`` gives the distribution over states after ``n`` steps starting from
    ``i``. ``P^0`` is the identity; every power is again row-stochastic.
    """
    size = _validate_matrix(P)
    if n < 0:
        raise ValueError("n must be non-negative")
    result = [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]
    base = [row[:] for row in P]
    e = n
    while e > 0:
        if e & 1:
            result = _matmul(result, base)
        base = _matmul(base, base)
        e >>= 1
    return result


def stationary_distribution(P, tol=1e-14, max_iter=100000):
    """Stationary distribution ``pi`` with ``pi P = pi`` (power iteration).

    Iterates a uniform start under ``P`` until convergence. For an irreducible
    aperiodic chain this is the unique long-run state distribution; the returned
    vector is non-negative and sums to one.
    """
    n = _validate_matrix(P)
    pi = [1.0 / n] * n
    for _ in range(max_iter):
        nxt = [sum(pi[i] * P[i][j] for i in range(n)) for j in range(n)]
        diff = sum(abs(nxt[j] - pi[j]) for j in range(n))
        pi = nxt
        if diff < tol:
            break
    s = sum(pi)
    return [x / s for x in pi]


def expected_hitting_time(P, target):
    """Expected number of steps to first reach ``target`` from each state.

    Solves the linear system ``h_i = 0`` for ``i = target`` and
    ``h_i = 1 + sum_j P_ij h_j`` otherwise, by Gaussian elimination. Returns the
    vector of expected hitting times (``0`` at the target, ``inf`` conceptually if
    unreachable -- the solver raises on a singular system in that case).
    """
    n = _validate_matrix(P)
    if not (0 <= target < n):
        raise ValueError("target out of range")
    # Build (I - Q) h = 1 for non-target states; h_target = 0.
    idx = [i for i in range(n) if i != target]
    m = len(idx)
    A = [[0.0] * m for _ in range(m)]
    b = [1.0] * m
    for r, i in enumerate(idx):
        for c, j in enumerate(idx):
            A[r][c] = (1.0 if i == j else 0.0) - P[i][j]
    # Gaussian elimination.
    for col in range(m):
        piv = max(range(col, m), key=lambda rr: abs(A[rr][col]))
        if abs(A[piv][col]) < 1e-15:
            raise ValueError("singular system (target may be unreachable)")
        A[col], A[piv] = A[piv], A[col]
        b[col], b[piv] = b[piv], b[col]
        pv = A[col][col]
        for j in range(col, m):
            A[col][j] /= pv
        b[col] /= pv
        for r in range(m):
            if r != col and A[r][col] != 0.0:
                f = A[r][col]
                for j in range(col, m):
                    A[r][j] -= f * A[col][j]
                b[r] -= f * b[col]
    h = [0.0] * n
    for r, i in enumerate(idx):
        h[i] = b[r]
    return h
