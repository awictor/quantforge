"""Sturm sequences: exact real-root counting and isolation for real polynomials.

A Sturm sequence of a polynomial ``p`` is ``p, p', `` then the negated remainders of
successive polynomial divisions. Sturm's theorem says the number of distinct real roots of
``p`` in ``(a, b]`` equals the drop in sign changes of the sequence evaluated at ``a`` versus
``b`` -- an exact integer count, no root-finding required. Bisecting on that count isolates
each real root in its own interval. Coefficients are highest-degree first, matching
:func:`quantforge.polyroots.polynomial_roots`. Pure standard library.
"""


def _deriv(p):
    n = len(p) - 1
    return [p[i] * (n - i) for i in range(n)] if n >= 1 else [0.0]


def _trim(p):
    i = 0
    while i < len(p) - 1 and abs(p[i]) < 1e-14:
        i += 1
    return p[i:]


def _poly_mod(a, b):
    # remainder of a divided by b (both highest-degree-first), float arithmetic
    r = _trim([float(x) for x in a])
    b = _trim(b)
    if len(b) == 1 and abs(b[0]) < 1e-14:
        return [0.0]
    # subtract shifted multiples of b until deg(r) < deg(b)
    while len(r) >= len(b) and not (len(r) == 1 and abs(r[0]) < 1e-14):
        factor = r[0] / b[0]
        shift = len(r) - len(b)
        for i in range(len(b)):
            r[i] -= factor * b[i]
        # r[0] is now ~0; drop it and continue
        r = _trim(r[1:] if len(r) > 1 else [0.0])
    return r


def sturm_sequence(coeffs):
    """Return the Sturm sequence of a real polynomial (list of coefficient lists).

    ``coeffs`` are highest-degree first. The sequence starts with the (square-free-agnostic)
    polynomial and its derivative; each subsequent term is the negated division remainder.
    """
    p = _trim([float(c) for c in coeffs])
    if len(p) == 1:
        return [p]
    seq = [p, _deriv(p)]
    while True:
        r = _poly_mod(seq[-2], seq[-1])
        r = _trim(r)
        if len(r) == 1 and abs(r[0]) < 1e-14:
            break
        seq.append([-c for c in r])
    return seq


def _eval(p, x):
    v = 0.0
    for c in p:
        v = v * x + c
    return v


def _sign_changes(seq, x):
    prev = 0
    changes = 0
    for p in seq:
        v = _eval(p, x)
        s = 0 if abs(v) < 1e-12 else (1 if v > 0 else -1)
        if s != 0:
            if prev != 0 and s != prev:
                changes += 1
            prev = s
    return changes


def real_root_count(coeffs, a, b):
    """Number of *distinct* real roots in the half-open interval ``(a, b]`` (Sturm's theorem).

    Exact for a square-free polynomial; repeated roots are counted once. ``a < b`` required.
    """
    if a >= b:
        raise ValueError("require a < b")
    seq = sturm_sequence(coeffs)
    return _sign_changes(seq, a) - _sign_changes(seq, b)


def isolate_real_roots(coeffs, a, b, max_depth=100):
    """Return disjoint sub-intervals of ``(a, b]``, each containing exactly one real root.

    Bisects the interval, using :func:`real_root_count` to decide which halves hold roots,
    until every returned interval isolates a single distinct real root. Each interval is a
    ``(lo, hi)`` pair.
    """
    if a >= b:
        raise ValueError("require a < b")
    seq = sturm_sequence(coeffs)

    def count(lo, hi):
        return _sign_changes(seq, lo) - _sign_changes(seq, hi)

    result = []
    stack = [(a, b, count(a, b), 0)]
    while stack:
        lo, hi, c, depth = stack.pop()
        if c == 0:
            continue
        if c == 1 or depth >= max_depth:
            result.append((lo, hi))
            continue
        mid = 0.5 * (lo + hi)
        cl = count(lo, mid)
        stack.append((mid, hi, c - cl, depth + 1))
        stack.append((lo, mid, cl, depth + 1))
    result.sort()
    return result
