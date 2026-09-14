"""Higher-order and derivative-free scalar root finders: Halley and secant.

`newton` in :mod:`quantforge.rootfind` converges quadratically; Halley's method adds the
second derivative to converge *cubically* (roughly tripling the correct digits per step),
and the secant method drops the derivative entirely, converging superlinearly (order ~1.618)
from two starting guesses. Both complement Newton where a second derivative is cheap
(Halley) or no derivative is available (secant). Pure standard library.
"""


def halley(f, fprime, fdoubleprime, x0, tol=1e-14, max_iter=100):
    """Halley's method: cubically-convergent root of ``f`` from ``x0``.

    Needs ``f``, its first derivative ``fprime``, and second derivative ``fdoubleprime``.
    The update is ``x - 2 f f' / (2 f'^2 - f f'')``. Raises if the denominator vanishes or
    it fails to converge in ``max_iter`` steps.
    """
    x = float(x0)
    for _ in range(max_iter):
        fx = f(x)
        if abs(fx) < tol:
            return x
        d1 = fprime(x)
        d2 = fdoubleprime(x)
        denom = 2 * d1 * d1 - fx * d2
        if denom == 0:
            raise ValueError("zero denominator in Halley step")
        x_next = x - 2 * fx * d1 / denom
        if abs(x_next - x) < tol:
            # accept only if it is genuinely a root, not a stagnated step
            if abs(f(x_next)) < 1e-6:
                return x_next
            raise ValueError("Halley's method stalled away from a root")
        x = x_next
    raise ValueError("Halley's method did not converge")


def secant(f, x0, x1, tol=1e-14, max_iter=200):
    """Secant method: derivative-free root of ``f`` from two initial guesses ``x0``, ``x1``.

    Approximates the derivative by the last two points; converges superlinearly
    (order ~1.618). Raises if the points coincide or it fails to converge.
    """
    x0 = float(x0)
    x1 = float(x1)
    f0 = f(x0)
    f1 = f(x1)
    for _ in range(max_iter):
        if abs(f1) < tol:
            return x1
        if f1 == f0:
            raise ValueError("secant step divides by zero (flat pair)")
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x2 - x1) < tol:
            if abs(f(x2)) < 1e-6:
                return x2
            raise ValueError("secant method stalled away from a root")
        x0, f0 = x1, f1
        x1, f1 = x2, f(x2)
    raise ValueError("secant method did not converge")
