"""Fast Walsh-Hadamard transform and bitwise (XOR/AND/OR) convolutions.

Where the FFT diagonalizes cyclic convolution, the Walsh-Hadamard transform diagonalizes
*XOR* convolution: ``(a * b)[k] = sum_{i ^ j = k} a[i] b[j]`` becomes a pointwise product in
the transform domain, computed in ``O(n log n)`` for a length that is a power of two.
Analogous transforms (the subset-sum / superset-sum zeta transforms) give AND and OR
convolutions. All are exact for integer inputs (no floating-point roots). Pure standard
library.
"""


def _check_pow2(a):
    n = len(a)
    if n == 0 or (n & (n - 1)) != 0:
        raise ValueError("length must be a positive power of two")
    return n


def fwht(a):
    """In-place-style fast Walsh-Hadamard transform (unnormalized). Returns a new list.

    Length must be a power of two. ``ifwht(fwht(a)) == a`` after the ``1/n`` scaling.
    """
    _check_pow2(a)
    r = list(a)
    n = len(r)
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                x, y = r[j], r[j + h]
                r[j] = x + y
                r[j + h] = x - y
        h *= 2
    return r


def ifwht(a):
    """Inverse Walsh-Hadamard transform (applies the ``1/n`` normalization)."""
    n = _check_pow2(a)
    r = fwht(a)
    return [v / n for v in r]


def xor_convolve(a, b):
    """XOR convolution: ``c[k] = sum_{i ^ j = k} a[i] b[j]``.

    ``a`` and ``b`` must share a power-of-two length. Integer inputs give exact integer
    output (the ``1/n`` division is exact because the Walsh transform's inverse sum is
    divisible by ``n`` for integer data).
    """
    n = _check_pow2(a)
    if len(b) != n:
        raise ValueError("inputs must have equal power-of-two length")
    fa = fwht(a)
    fb = fwht(b)
    fc = [fa[i] * fb[i] for i in range(n)]
    conv = fwht(fc)
    out = []
    for v in conv:
        q, rem = divmod(v, n)
        out.append(q if rem == 0 else v / n)
    return out


def _zeta_subset(a):
    # sum over subsets: f[mask] = sum_{sub subset of mask} a[sub]
    n = _check_pow2(a)
    r = list(a)
    bit = 1
    while bit < n:
        for mask in range(n):
            if mask & bit:
                r[mask] += r[mask ^ bit]
        bit <<= 1
    return r


def _zeta_superset(a):
    # sum over supersets: f[mask] = sum_{sup superset of mask} a[sup]
    n = _check_pow2(a)
    r = list(a)
    bit = 1
    while bit < n:
        for mask in range(n):
            if not (mask & bit):
                r[mask] += r[mask | bit]
        bit <<= 1
    return r


def and_convolve(a, b):
    """AND convolution: ``c[k] = sum_{i & j = k} a[i] b[j]`` (via superset zeta transform)."""
    n = _check_pow2(a)
    if len(b) != n:
        raise ValueError("inputs must have equal power-of-two length")
    fa = _zeta_superset(a)
    fb = _zeta_superset(b)
    fc = [fa[i] * fb[i] for i in range(n)]
    # inverse superset zeta (Mobius): subtract supersets back out
    r = list(fc)
    bit = 1
    while bit < n:
        for mask in range(n):
            if not (mask & bit):
                r[mask] -= r[mask | bit]
        bit <<= 1
    return r


def or_convolve(a, b):
    """OR convolution: ``c[k] = sum_{i | j = k} a[i] b[j]`` (via subset zeta transform)."""
    n = _check_pow2(a)
    if len(b) != n:
        raise ValueError("inputs must have equal power-of-two length")
    fa = _zeta_subset(a)
    fb = _zeta_subset(b)
    fc = [fa[i] * fb[i] for i in range(n)]
    # inverse subset zeta (Mobius): subtract subsets back out
    r = list(fc)
    bit = 1
    while bit < n:
        for mask in range(n):
            if mask & bit:
                r[mask] -= r[mask ^ bit]
        bit <<= 1
    return r
