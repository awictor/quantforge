"""One-dimensional wave equation solver (explicit finite difference).

Solves ``u_tt = c^2 u_xx`` on ``[0, L]`` with fixed (Dirichlet) ends, an initial shape
``u(x, 0)`` and initial velocity ``u_t(x, 0)``. The explicit central-difference
(leapfrog) scheme steps the profile forward using two time levels; it is second-order
accurate and stable while the Courant number ``C = c dt / dx <= 1`` (the CFL condition).
Standing waves oscillate, and a localized pulse splits into left- and right-moving
halves per d'Alembert. Pure standard library.
"""


def wave_equation(u0, v0, c, dx, dt, n_steps, left=0.0, right=0.0):
    """Explicit finite-difference evolution of ``u_tt = c^2 u_xx``.

    ``u0`` is the initial displacement profile (grid values), ``v0`` the initial
    velocity (same length; pass zeros for a plucked start). Fixed ends at ``left`` /
    ``right``. Returns the displacement after ``n_steps`` steps. Requires the CFL
    condition ``c dt / dx <= 1`` for stability (raises otherwise).
    """
    n = len(u0)
    if n < 3:
        raise ValueError("need at least 3 grid points")
    if len(v0) != n:
        raise ValueError("u0 and v0 must have equal length")
    if c <= 0 or dx <= 0 or dt <= 0:
        raise ValueError("c, dx, dt must be positive")
    C2 = (c * dt / dx) ** 2
    if C2 > 1.0 + 1e-12:
        raise ValueError("CFL violated: need c*dt/dx <= 1")

    prev = list(u0)
    prev[0], prev[-1] = left, right
    # First step uses the initial velocity (Taylor): u^1 = u^0 + dt v0 + 0.5 C2 (u_xx).
    cur = [0.0] * n
    cur[0], cur[-1] = left, right
    for i in range(1, n - 1):
        cur[i] = (prev[i] + dt * v0[i]
                  + 0.5 * C2 * (prev[i + 1] - 2 * prev[i] + prev[i - 1]))
    if n_steps == 1:
        return cur

    for _ in range(n_steps - 1):
        nxt = [0.0] * n
        nxt[0], nxt[-1] = left, right
        for i in range(1, n - 1):
            nxt[i] = (2 * cur[i] - prev[i]
                      + C2 * (cur[i + 1] - 2 * cur[i] + cur[i - 1]))
        prev, cur = cur, nxt
    return cur
