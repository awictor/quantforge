"""Remez exchange algorithm: best (minimax) polynomial approximation.

Least-squares fits (:func:`quantforge.chebyshev.chebyshev_fit`) minimize the *average* squared
error. The Remez algorithm instead minimizes the *worst-case* error ``max_x |f(x) - p(x)|`` --
the true minimax polynomial. By Chebyshev's equioscillation theorem the optimal degree-``n``
approximation has an error curve that alternates sign and reaches its maximum magnitude at
``n + 2`` points. Remez finds it by iteration: solve for the polynomial that makes the error
equal-and-alternating at the current reference points, then move each reference to the local
extremum of the resulting error, and repeat until the extrema equalize.

Returns the polynomial coefficients (highest degree first) and the achieved minimax error. Pure
standard library.
"""

from .lu import lu_solve


def _poly_eval(coeffs, x):
    # coeffs highest-degree first
    r = 0.0
    for c in coeffs:
        r = r * x + c
    return r


def remez(f, a, b, degree, iterations=100, tol=1e-12, grid=2000):
    """Minimax degree-``degree`` polynomial approximation of ``f`` on ``[a, b]``.

    Returns a dict with ``coeffs`` (highest-degree first, length ``degree + 1``), ``error`` (the
    minimax error ``max |f - p|``), and ``n_iter``. Uses the Remez exchange with a dense grid to
    locate error extrema.
    """
    n = degree
    m = n + 2                                # number of reference points
    # initial references: Chebyshev extrema mapped to [a, b]
    import math
    refs = [0.5 * (a + b) + 0.5 * (b - a) * math.cos(math.pi * i / (m - 1)) for i in range(m)]
    refs.sort()

    coeffs = [0.0] * (n + 1)
    err = 0.0
    for it in range(iterations):
        # Solve the linear system: for each reference x_i,
        #   p(x_i) + (-1)^i * E = f(x_i),  unknowns = coeffs (n+1) and E.
        A = []
        rhs = []
        for i, x in enumerate(refs):
            row = [x ** (n - k) for k in range(n + 1)] + [(-1.0) ** i]
            A.append(row)
            rhs.append(f(x))
        sol = lu_solve(A, rhs)
        coeffs = sol[:n + 1]
        E = sol[n + 1]

        # Find the grid point of maximum error and the local extrema for the exchange.
        best_x = a
        best_abs = -1.0
        errs = []
        for j in range(grid + 1):
            x = a + (b - a) * j / grid
            e = f(x) - _poly_eval(coeffs, x)
            errs.append((x, e))
            if abs(e) > best_abs:
                best_abs = abs(e)
                best_x = x

        # Rebuild the reference set from the alternating local extrema of the error curve.
        extrema = _alternating_extrema(errs, m)
        if extrema:
            refs = extrema
        # convergence: max error close to |E|
        if abs(best_abs - abs(E)) < tol * max(1.0, best_abs):
            err = best_abs
            break
        err = best_abs
    return {"coeffs": coeffs, "error": err, "n_iter": it + 1}


def _alternating_extrema(errs, m):
    # collect local extrema (sign-consistent runs), pick m with alternating signs, largest each
    extrema = []
    for i in range(1, len(errs) - 1):
        x, e = errs[i]
        el = errs[i - 1][1]
        er = errs[i + 1][1]
        if (e >= el and e >= er) or (e <= el and e <= er):
            extrema.append((x, e))
    # include endpoints
    extrema = [errs[0]] + extrema + [errs[-1]]
    # merge consecutive same-sign extrema, keeping the larger magnitude
    merged = []
    for x, e in extrema:
        if merged and (e > 0) == (merged[-1][1] > 0):
            if abs(e) > abs(merged[-1][1]):
                merged[-1] = (x, e)
        else:
            merged.append((x, e))
    if len(merged) < m:
        return None
    # keep the m largest-magnitude, then sort by x
    merged.sort(key=lambda p: -abs(p[1]))
    keep = sorted(p[0] for p in merged[:m])
    return keep
