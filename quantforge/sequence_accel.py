"""Sequence acceleration: Aitken's delta-squared, Shanks transform, Steffensen.

A linearly-converging sequence (fixed-point iterations, slowly-summing series) can be
sped up dramatically by extrapolating its limit from the pattern of successive terms:

  * ``aitken`` -- Aitken's ``Delta^2`` process on a whole sequence, returning the
    accelerated sequence ``s* = s_n - (s_{n+1} - s_n)^2 / (s_{n+2} - 2 s_{n+1} + s_n)``.
  * ``shanks`` -- the Shanks transform, algebraically the same extrapolation stated
    per triple; iterating it repeatedly is Aitken applied again and again.
  * ``steffensen`` -- Steffensen's method: Aitken acceleration fused into a fixed-point
    iteration ``x = g(x)``, giving quadratic convergence from a linear map without any
    derivative.

Pure standard library.
"""


def _accel_point(s0, s1, s2):
    """One Aitken/Shanks extrapolation from three consecutive terms."""
    denom = s2 - 2.0 * s1 + s0
    if denom == 0.0:
        return s2
    return s0 - (s1 - s0) ** 2 / denom


def aitken(sequence):
    """Aitken's delta-squared acceleration of a sequence.

    Returns the accelerated sequence, two elements shorter than the input. Each output
    term extrapolates the limit from three consecutive input terms; for a linearly
    convergent sequence it converges markedly faster.
    """
    s = list(sequence)
    if len(s) < 3:
        raise ValueError("need at least 3 terms")
    return [_accel_point(s[i], s[i + 1], s[i + 2]) for i in range(len(s) - 2)]


def shanks(sequence):
    """Shanks transform of a sequence (one pass; same formula as :func:`aitken`)."""
    return aitken(sequence)


def steffensen(g, x0, tol=1e-12, max_iter=100):
    """Steffensen's method: quadratically-convergent fixed point of ``x = g(x)``.

    Applies Aitken acceleration to the fixed-point iterates, so it converges
    quadratically like Newton's method but needs no derivative -- only the map ``g``.
    Returns a dict with the ``root`` (the fixed point), ``iterations`` and
    ``converged``. Raises if a zero denominator stalls the iteration before
    convergence.
    """
    x = float(x0)
    for it in range(1, max_iter + 1):
        gx = g(x)
        ggx = g(gx)
        denom = ggx - 2.0 * gx + x
        if denom == 0.0:
            if abs(gx - x) <= tol:
                return {"root": gx, "iterations": it, "converged": True}
            raise ValueError("zero denominator in Steffensen iteration")
        x_new = x - (gx - x) ** 2 / denom
        if abs(x_new - x) <= tol:
            return {"root": x_new, "iterations": it, "converged": True}
        x = x_new
    return {"root": x, "iterations": max_iter, "converged": False}
