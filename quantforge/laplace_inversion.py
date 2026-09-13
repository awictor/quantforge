"""Numerical inversion of a Laplace transform (Gaver-Stehfest).

Given a Laplace transform ``F(s)`` known only as a function (e.g. a closed-form
transform of an unknown time function), recover the time-domain ``f(t)``. The
Gaver-Stehfest algorithm evaluates ``F`` at a set of real points ``k ln2 / t`` and
combines them with tabulated coefficients,

    f(t) ~ (ln2 / t) sum_{k=1}^{N} V_k F(k ln2 / t),

using only real arithmetic (no complex contour). It is exact-ish for smooth,
non-oscillatory ``f`` -- the workhorse for inverting transforms in queueing, finance
(e.g. arithmetic-Asian and occupation-time formulas), and diffusion problems. ``N`` even
(typically 8-16); higher ``N`` needs higher working precision, so accuracy peaks and
then degrades. Pure standard library.
"""

import math


def _stehfest_coefficients(N):
    """Salzer coefficients V_k for the Gaver-Stehfest method (N even)."""
    if N % 2 != 0:
        raise ValueError("N must be even")
    half = N // 2
    V = [0.0] * (N + 1)          # 1-based
    for k in range(1, N + 1):
        s = 0.0
        jlo = (k + 1) // 2
        jhi = min(k, half)
        for j in range(jlo, jhi + 1):
            num = (j ** half) * math.factorial(2 * j)
            den = (math.factorial(half - j) * math.factorial(j)
                   * math.factorial(j - 1) * math.factorial(k - j)
                   * math.factorial(2 * j - k))
            s += num / den
        V[k] = ((-1) ** (k + half)) * s
    return V


def laplace_inversion(F, t, N=12):
    """Invert the Laplace transform ``F(s)`` at time ``t`` by Gaver-Stehfest.

    ``F`` is a callable of a single real argument ``s > 0``. ``N`` (even, default 12) is
    the number of terms. Returns the estimated ``f(t)``. Best for smooth,
    non-oscillatory functions; oscillatory or discontinuous ``f`` invert poorly.
    """
    if t <= 0:
        raise ValueError("t must be positive")
    if N % 2 != 0:
        raise ValueError("N must be even")
    V = _stehfest_coefficients(N)
    ln2_over_t = math.log(2.0) / t
    total = 0.0
    for k in range(1, N + 1):
        total += V[k] * F(k * ln2_over_t)
    return ln2_over_t * total
