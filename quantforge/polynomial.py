"""Dense polynomial arithmetic on coefficient lists (low-degree-first).

A small, exact polynomial algebra: add, subtract, multiply, divide with remainder,
differentiate, integrate, evaluate (Horner), and the polynomial GCD. Coefficients are
ordered low-degree-first -- ``[c0, c1, c2]`` means ``c0 + c1 x + c2 x^2`` -- the
natural order for calculus operations. These complement :func:`quantforge.polynomial_roots`
(all roots) and the FFT-based :func:`quantforge.convolve` (fast product for long
polynomials). Pure standard library.
"""


def _trim(c):
    """Drop trailing (high-degree) near-zero coefficients; keep at least [0]."""
    out = list(c)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out if out else [0.0]


def poly_add(a, b):
    """Sum of two polynomials (coefficient lists, low-degree-first)."""
    n = max(len(a), len(b))
    return _trim([(a[i] if i < len(a) else 0.0) + (b[i] if i < len(b) else 0.0)
                  for i in range(n)])


def poly_sub(a, b):
    """Difference ``a - b`` of two polynomials."""
    n = max(len(a), len(b))
    return _trim([(a[i] if i < len(a) else 0.0) - (b[i] if i < len(b) else 0.0)
                  for i in range(n)])


def poly_mul(a, b):
    """Product of two polynomials by direct convolution.

    Exact for the small polynomials typical of algebra; for long polynomials the
    FFT-based :func:`quantforge.convolve` is asymptotically faster.
    """
    if not a or not b:
        return [0.0]
    out = [0.0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return _trim(out)


def poly_divmod(num, den):
    """Polynomial long division: return ``(quotient, remainder)``.

    ``num = quotient * den + remainder`` with ``deg(remainder) < deg(den)``. Raises on
    a zero divisor.
    """
    d = _trim(den)
    if d == [0.0]:
        raise ValueError("division by the zero polynomial")
    r = _trim(num)
    dd = len(d) - 1
    lead = d[-1]
    quo = [0.0] * max(len(r) - dd, 1)
    while len(r) - 1 >= dd and r != [0.0]:
        shift = len(r) - 1 - dd
        factor = r[-1] / lead
        quo[shift] = factor
        # Subtract factor * x^shift * d from r. The top coefficient cancels exactly
        # by construction, but float round-off can leave a tiny residue that _trim's
        # exact == 0 test would miss (causing a non-terminating loop), so force it.
        for i, di in enumerate(d):
            r[shift + i] -= factor * di
        r[-1] = 0.0
        r = _trim(r)
    return _trim(quo), r


def poly_derivative(c):
    """Derivative of a polynomial: ``[c1, 2 c2, 3 c3, ...]``."""
    if len(c) <= 1:
        return [0.0]
    return [i * c[i] for i in range(1, len(c))]


def poly_integral(c, constant=0.0):
    """Antiderivative of a polynomial, with integration constant ``constant``."""
    return [constant] + [c[i] / (i + 1) for i in range(len(c))]


def poly_eval(c, x):
    """Evaluate a polynomial at ``x`` by Horner's method."""
    r = 0.0
    for coef in reversed(c):
        r = r * x + coef
    return r


def poly_gcd(a, b, tol=1e-9):
    """Monic greatest common divisor of two polynomials (Euclidean algorithm).

    Returns the GCD normalized to a monic polynomial (leading coefficient 1); useful
    for detecting and factoring out repeated roots (``gcd(p, p')``). Coefficients below
    ``tol`` in the remainder are treated as zero to tame round-off.
    """
    x = _trim(a)
    y = _trim(b)
    while y != [0.0]:
        _, r = poly_divmod(x, y)
        r = _trim([0.0 if abs(v) < tol else v for v in r])
        x, y = y, r
    # Normalize to monic.
    lead = x[-1]
    if lead != 0:
        x = [v / lead for v in x]
    return _trim(x)
