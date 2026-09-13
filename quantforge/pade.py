"""Pade approximants and continued-fraction evaluation.

A Pade approximant ``[m/n]`` is the ratio of a degree-``m`` over a degree-``n``
polynomial whose Taylor expansion matches a given series through order ``m + n``. It
often models functions -- especially near singularities or over wide ranges -- far
better than the truncated Taylor series it is built from, converging where the series
diverges. This module builds the approximant from Taylor coefficients and evaluates it,
and provides a general modified-Lentz continued-fraction evaluator.

Pure standard library.
"""


def pade(coeffs, m, n):
    """Pade ``[m/n]`` approximant from Taylor coefficients ``coeffs``.

    ``coeffs[k]`` is the coefficient of ``x^k``; needs at least ``m + n + 1`` of them.
    Returns ``(num, den)`` -- the numerator (length ``m+1``) and denominator (length
    ``n+1``, normalized to ``den[0] = 1``) coefficient lists. The approximant's Taylor
    series matches ``coeffs`` through order ``m + n``.
    """
    if m < 0 or n < 0:
        raise ValueError("m and n must be non-negative")
    if len(coeffs) < m + n + 1:
        raise ValueError("need at least m + n + 1 Taylor coefficients")
    c = [float(coeffs[k]) for k in range(m + n + 1)]

    if n == 0:
        return c[:m + 1], [1.0]

    # Solve for denominator coefficients b_1..b_n from the n equations
    # sum_{j=1..n} c[m+i-j] b_j = -c[m+i], i = 1..n.
    A = [[c[m + i - j] if 0 <= m + i - j else 0.0 for j in range(1, n + 1)]
         for i in range(1, n + 1)]
    rhs = [-c[m + i] for i in range(1, n + 1)]
    b = _solve(A, rhs)
    den = [1.0] + b

    # Numerator a_k = sum_{j=0..min(k,n)} b_j c[k-j], with b_0 = 1.
    num = []
    for k in range(m + 1):
        s = c[k]
        for j in range(1, min(k, n) + 1):
            s += b[j - 1] * c[k - j]
        num.append(s)
    return num, den


def _solve(A, rhs):
    """Small dense linear solve by Gaussian elimination with partial pivoting."""
    n = len(A)
    M = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-300:
            raise ValueError("singular Pade system (try different m, n)")
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        for r in range(n):
            if r != col:
                f = M[r][col] / pv
                for k in range(col, n + 1):
                    M[r][k] -= f * M[col][k]
    return [M[i][n] / M[i][i] for i in range(n)]


def pade_eval(num, den, x):
    """Evaluate a Pade approximant ``(num, den)`` at ``x`` by Horner's method."""
    p = 0.0
    for a in reversed(num):
        p = p * x + a
    q = 0.0
    for b in reversed(den):
        q = q * x + b
    if q == 0.0:
        raise ValueError("Pade denominator vanishes at x")
    return p / q


def lentz_continued_fraction(a, b, tol=1e-15, max_iter=1000, tiny=1e-300):
    """Evaluate a continued fraction by the modified Lentz algorithm.

    Computes ``b0 + a1/(b1 + a2/(b2 + ...))`` where ``a(k)`` and ``b(k)`` are callables
    giving the ``k``-th partial numerator and denominator (``b(0)`` is the leading
    term, ``a(0)`` is unused). Returns the converged value. Robust to zero
    intermediate values via the ``tiny`` guard.
    """
    f = b(0)
    if f == 0.0:
        f = tiny
    c = f
    d = 0.0
    for k in range(1, max_iter + 1):
        ak = a(k)
        bk = b(k)
        d = bk + ak * d
        if d == 0.0:
            d = tiny
        c = bk + ak / c
        if c == 0.0:
            c = tiny
        d = 1.0 / d
        delta = c * d
        f *= delta
        if abs(delta - 1.0) < tol:
            return f
    return f
