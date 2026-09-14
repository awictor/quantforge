"""Polynomial arithmetic and Lagrange interpolation over a prime field.

Working with polynomials modulo a prime ``p`` -- exact, no floating error -- underpins
secret sharing, error-correcting codes, and the fast-transform machinery. This module gives
evaluation, addition, and multiplication of coefficient lists mod ``p``, and Lagrange
interpolation to recover the unique degree-``< n`` polynomial through ``n`` points
``(x_i, y_i)`` with distinct ``x_i``. Coefficients are lowest-degree first here (``c[0]`` is
the constant term). Pure standard library.
"""


def poly_eval_mod(coeffs, x, mod):
    """Evaluate a polynomial (lowest-degree-first) at ``x`` modulo ``mod`` (Horner)."""
    result = 0
    for c in reversed(coeffs):
        result = (result * x + c) % mod
    return result


def poly_add_mod(a, b, mod):
    """Sum of two polynomials (lowest-degree-first) modulo ``mod``."""
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(len(a)):
        out[i] = a[i] % mod
    for i in range(len(b)):
        out[i] = (out[i] + b[i]) % mod
    return _trim_mod(out)


def poly_mul_mod(a, b, mod):
    """Product of two polynomials (lowest-degree-first) modulo ``mod``."""
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] = (out[i + j] + ai * bj) % mod
    return _trim_mod(out)


def _trim_mod(c):
    c = list(c)
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


def lagrange_interpolate_mod(points, mod):
    """Recover the polynomial through ``points`` = ``[(x_i, y_i)]`` modulo prime ``mod``.

    Returns coefficients (lowest-degree first) of the unique polynomial of degree ``< n``
    that passes through all ``n`` points; ``x_i`` must be distinct modulo ``mod``. Builds
    ``sum_i y_i * prod_{j!=i} (x - x_j)/(x_i - x_j)`` with exact modular inverses.
    """
    xs = [x % mod for x, _ in points]
    ys = [y % mod for _, y in points]
    n = len(points)
    if len(set(xs)) != n:
        raise ValueError("x-coordinates must be distinct modulo mod")
    result = [0]
    for i in range(n):
        # numerator polynomial prod_{j != i} (x - x_j)
        num = [1]
        denom = 1
        for j in range(n):
            if j == i:
                continue
            num = poly_mul_mod(num, [(-xs[j]) % mod, 1], mod)
            denom = (denom * (xs[i] - xs[j])) % mod
        inv = pow(denom, mod - 2, mod)
        scale = (ys[i] * inv) % mod
        term = [(c * scale) % mod for c in num]
        result = poly_add_mod(result, term, mod)
    return result if result != [0] else [0]
