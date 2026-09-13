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


def cumulative_default_term_structure(P, default_state, horizons,
                                      start_state=0):
    """Cumulative default probability by horizon from a rating-migration matrix.

    ``P`` is a one-period rating transition matrix with ``default_state`` an
    absorbing default row. For each ``n`` in ``horizons`` the cumulative default
    probability from ``start_state`` is the default-column entry of ``P^n``, i.e.
    ``(P^n)[start_state][default_state]``. Non-decreasing in the horizon (default
    is absorbing) and rising toward one if default is reachable.
    """
    n_states = _validate_matrix(P)
    if not (0 <= default_state < n_states and 0 <= start_state < n_states):
        raise ValueError("state indices out of range")
    if abs(P[default_state][default_state] - 1.0) > 1e-9:
        raise ValueError("default_state must be absorbing (self-transition 1)")
    out = []
    for n in horizons:
        Pn = n_step_transition(P, n)
        out.append(Pn[start_state][default_state])
    return out


def marginal_default_probabilities(P, default_state, horizons, start_state=0):
    """Marginal (per-period) default probabilities between successive horizons.

    Differences of the :func:`cumulative_default_term_structure`; each is the
    probability of defaulting in ``(horizons[k-1], horizons[k]]`` having survived
    to ``horizons[k-1]``. Non-negative because the cumulative curve is
    non-decreasing.
    """
    cum = cumulative_default_term_structure(P, default_state, horizons,
                                            start_state)
    out = []
    prev = 0.0
    for c in cum:
        out.append(c - prev)
        prev = c
    return out


def _invert(M):
    """Invert a small dense matrix by Gauss-Jordan elimination."""
    n = len(M)
    a = [row[:] + [1.0 if i == j else 0.0 for j in range(n)]
         for i, row in enumerate(M)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) < 1e-15:
            raise ValueError("singular matrix")
        a[col], a[piv] = a[piv], a[col]
        pv = a[col][col]
        for j in range(2 * n):
            a[col][j] /= pv
        for r in range(n):
            if r != col:
                f = a[r][col]
                for j in range(2 * n):
                    a[r][j] -= f * a[col][j]
    return [row[n:] for row in a]


def fundamental_matrix(P, transient_states):
    """Fundamental matrix ``N = (I - Q)^{-1}`` of an absorbing chain.

    ``transient_states`` lists the indices of the non-absorbing states; ``Q`` is
    their sub-transition block. ``N_ij`` is the expected number of visits to
    transient state ``j`` starting from ``i`` before absorption. Requires the chain
    to be absorbing (every transient state eventually reaches an absorbing one).
    """
    _validate_matrix(P)
    m = len(transient_states)
    if m == 0:
        raise ValueError("need at least one transient state")
    Q = [[P[i][j] for j in transient_states] for i in transient_states]
    I_minus_Q = [[(1.0 if r == c else 0.0) - Q[r][c] for c in range(m)]
                 for r in range(m)]
    return _invert(I_minus_Q)


def expected_steps_to_absorption(P, transient_states):
    """Expected steps to absorption from each transient state.

    Row sums of the :func:`fundamental_matrix` ``N`` -- the total expected visits
    across all transient states before hitting an absorbing state. Positive for
    every transient state.
    """
    N = fundamental_matrix(P, transient_states)
    return [sum(row) for row in N]


def absorption_probabilities(P, transient_states, absorbing_states):
    """Probability of ending in each absorbing state from each transient state.

    ``B = N R`` where ``N`` is the :func:`fundamental_matrix` and ``R`` is the
    transient-to-absorbing transition block. Row ``i`` (a transient state) is a
    distribution over ``absorbing_states`` summing to one.
    """
    N = fundamental_matrix(P, transient_states)
    m = len(transient_states)
    k = len(absorbing_states)
    R = [[P[transient_states[r]][absorbing_states[c]] for c in range(k)]
         for r in range(m)]
    return [[sum(N[i][t] * R[t][c] for t in range(m)) for c in range(k)]
            for i in range(m)]


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


def _validate_generator(Q):
    """Check ``Q`` is a valid continuous-time generator: rows sum ~0, off-diag >= 0."""
    n = len(Q)
    if n == 0 or any(len(row) != n for row in Q):
        raise ValueError("Q must be a non-empty square matrix")
    for i in range(n):
        if abs(sum(Q[i])) > 1e-6:
            raise ValueError("generator rows must sum to zero")
        for j in range(n):
            if i != j and Q[i][j] < -1e-12:
                raise ValueError("off-diagonal generator entries must be non-negative")
        if Q[i][i] > 1e-12:
            raise ValueError("diagonal generator entries must be non-positive")
    return n


def generator_to_transition(Q, t=1.0):
    """Transition matrix ``P(t) = exp(Q t)`` of a continuous-time Markov chain.

    ``Q`` is a rate generator (rows summing to zero, non-negative off-diagonals).
    Returns the ``t``-horizon transition matrix, whose rows sum to one with
    non-negative entries. ``P(0)`` is the identity and ``P(s) P(t) = P(s + t)``
    (the semigroup property). The basis for continuous-time rating migration.
    """
    from .matrix_exp import matrix_exp
    _validate_generator(Q)
    if t < 0.0:
        raise ValueError("t must be non-negative")
    n = len(Q)
    Qt = [[Q[i][j] * t for j in range(n)] for i in range(n)]
    P = matrix_exp(Qt)
    # Clean tiny negatives and renormalize rows against roundoff.
    for i in range(n):
        row = [max(0.0, P[i][j]) for j in range(n)]
        s = sum(row)
        P[i] = [v / s for v in row] if s > 0 else row
    return P


def generator_default_probability(Q, default_state, horizons, start_state=0):
    """Cumulative default probability at each horizon from a rating generator.

    Exponentiates ``Q`` to each horizon and reads the transition probability from
    ``start_state`` into the absorbing ``default_state``. Returns one probability per
    horizon; non-decreasing when the default state is absorbing.
    """
    _validate_generator(Q)
    out = []
    for h in horizons:
        P = generator_to_transition(Q, h)
        out.append(P[start_state][default_state])
    return out
