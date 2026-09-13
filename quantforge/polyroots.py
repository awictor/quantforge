"""All roots of a polynomial by the Durand-Kerner (Weierstrass) iteration.

Given real or complex coefficients, this finds *every* root at once -- real and
complex -- without deflation. Each root estimate is refined by

    z_i <- z_i - p(z_i) / prod_{j != i} (z_i - z_j),

a simultaneous Newton-like step that converges from spread-out complex starting points
(the classic choice ``(0.4 + 0.9i)^k``). No derivative of ``p`` beyond ``p`` itself is
needed, and it handles the full complex spectrum a real bisection/Newton root-finder
cannot. Pure standard library.
"""


def _horner(coeffs, z):
    """Evaluate a polynomial (coeffs high-degree first) at ``z``."""
    r = 0j
    for c in coeffs:
        r = r * z + c
    return r


def polynomial_roots(coeffs, tol=1e-12, max_iter=500):
    """All roots of a polynomial by the Durand-Kerner method.

    ``coeffs`` are the coefficients from the highest degree down (e.g. ``[1, -3, 2]``
    for ``x^2 - 3x + 2``); real or complex. Returns a list of the ``n`` roots as
    complex numbers (a root with a negligible imaginary part is still returned as
    ``complex`` -- take ``.real`` if you know it is real). Leading zeros are trimmed.
    """
    # Trim leading zeros.
    c = list(coeffs)
    while len(c) > 1 and c[0] == 0:
        c.pop(0)
    n = len(c) - 1
    if n < 1:
        raise ValueError("need a polynomial of degree >= 1")

    # Normalize to monic.
    lead = c[0]
    c = [ci / lead for ci in c]

    # Spread-out complex starting points: (0.4 + 0.9i)^k.
    seed = complex(0.4, 0.9)
    roots = [seed ** k for k in range(n)]

    for _ in range(max_iter):
        max_delta = 0.0
        new = list(roots)
        for i in range(n):
            zi = roots[i]
            denom = 1.0 + 0j
            for j in range(n):
                if j != i:
                    denom *= (zi - roots[j])
            if denom == 0:
                continue
            delta = _horner(c, zi) / denom
            new[i] = zi - delta
            d = abs(delta)
            if d > max_delta:
                max_delta = d
        roots = new
        if max_delta < tol:
            break
    return roots
