"""One-dimensional heat/diffusion equation solver (Crank-Nicolson).

Solves ``u_t = alpha u_xx`` on ``[0, L]`` with Dirichlet boundary values and a given
initial profile, by the Crank-Nicolson scheme -- the average of explicit and implicit
Euler in time. It is second-order accurate in both space and time and unconditionally
stable (no CFL step-size restriction), so it handles diffusive smoothing, heat
conduction, and the transformed Black-Scholes equation robustly. Each step solves a
tridiagonal system. Pure standard library.
"""

from .andreasenhuge import _solve_tridiag


def heat_equation_cn(u0, alpha, dx, dt, n_steps, left=None, right=None):
    """Crank-Nicolson evolution of ``u_t = alpha u_xx``.

    ``u0`` is the initial profile (interior + boundary grid values). ``dx``/``dt`` the
    space/time steps, ``alpha`` the diffusivity. ``left``/``right`` fix the Dirichlet
    boundary values (default: hold the initial endpoints). Returns the profile after
    ``n_steps`` steps. Unconditionally stable, second-order in space and time.
    """
    n = len(u0)
    if n < 3:
        raise ValueError("need at least 3 grid points")
    if alpha < 0 or dx <= 0 or dt <= 0:
        raise ValueError("alpha, dx, dt must be positive")
    u = list(u0)
    bl = u[0] if left is None else left
    br = u[-1] if right is None else right
    u[0], u[-1] = bl, br

    r = alpha * dt / (2.0 * dx * dx)
    m = n - 2                      # interior unknowns
    # Implicit-side tridiagonal (I - r L): diag 1+2r, off -r.
    lower = [-r] * m
    diag = [1.0 + 2.0 * r] * m
    upper = [-r] * m

    for _ in range(n_steps):
        # RHS = (I + r L) u^k on the interior, plus boundary contributions.
        rhs = []
        for i in range(1, n - 1):
            val = (1.0 - 2.0 * r) * u[i] + r * (u[i - 1] + u[i + 1])
            rhs.append(val)
        # Boundary terms from the implicit side (constant Dirichlet values).
        rhs[0] += r * bl
        rhs[-1] += r * br
        interior = _solve_tridiag(list(lower), list(diag), list(upper), rhs)
        u = [bl] + interior + [br]
    return u
