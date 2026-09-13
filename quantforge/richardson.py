"""General Richardson extrapolation to the zero-step-size limit.

Many numerical estimates ``A(h)`` converge to the true value as ``A(h) = A + c h^p +
O(h^q)``. Richardson extrapolation combines estimates at successively halved steps to
cancel the leading error term, then the next, building a triangular tableau whose
diagonal converges far faster than the raw sequence. Given estimates ``A(h), A(h/2),
A(h/4), ...`` and the leading error order ``p`` (and step ratio ``t``, usually 2), it
returns the accelerated limit -- the engine behind Romberg integration and Richardson
differentiation, here as a standalone tool for any such sequence. Pure standard library.
"""


def richardson_extrapolate(estimates, p=1.0, t=2.0):
    """Richardson-extrapolate a sequence of step-halved estimates to the ``h -> 0`` limit.

    ``estimates[i] = A(h / t^i)``, ordered from coarsest to finest. ``p`` is the leading
    error exponent (1 for first-order, 2 for a central difference, etc.), ``t`` the step
    ratio between consecutive estimates. Returns the best (last-diagonal) extrapolated
    value. Successive columns cancel the ``h^p, h^{p+1}, ...`` error terms.
    """
    n = len(estimates)
    if n == 0:
        raise ValueError("need at least one estimate")
    # Neville-style tableau: T[i][j] cancels j leading error terms.
    T = [[0.0] * n for _ in range(n)]
    for i in range(n):
        T[i][0] = estimates[i]
    for j in range(1, n):
        factor = t ** (p * j)
        for i in range(j, n):
            T[i][j] = (factor * T[i][j - 1] - T[i - 1][j - 1]) / (factor - 1.0)
    return T[n - 1][n - 1]


def richardson_table(estimates, p=1.0, t=2.0):
    """Full Richardson tableau (list of rows) for inspecting convergence.

    Row ``i`` holds ``T[i][0..i]``; the diagonal ``T[i][i]`` is the order-``i``
    extrapolation. Handy to watch the estimate stabilize down the diagonal.
    """
    n = len(estimates)
    if n == 0:
        raise ValueError("need at least one estimate")
    T = [[estimates[i]] for i in range(n)]
    for i in range(1, n):
        for j in range(1, i + 1):
            factor = t ** (p * j)
            T[i].append((factor * T[i][j - 1] - T[i - 1][j - 1]) / (factor - 1.0))
    return T
