"""Symplectic integrators for Hamiltonian systems (velocity Verlet, leapfrog).

For a system ``q' = p/m``, ``p' = F(q)`` (a separable Hamiltonian), a general RK method
slowly drifts the total energy over long integrations. Symplectic integrators preserve
the geometric structure of phase space, so the energy stays *bounded* -- oscillating
about the true value rather than drifting -- making them the right tool for orbital
mechanics, molecular dynamics, and any long-horizon conservative simulation. Velocity
Verlet is second-order and time-reversible. ``q`` and ``p`` may be scalars or vectors.
Pure standard library.
"""


def _vec(x):
    return list(x) if isinstance(x, (list, tuple)) else [x]


def velocity_verlet(force, q0, p0, mass, dt, n_steps):
    """Velocity-Verlet integration of ``q'' = force(q) / mass``.

    ``force(q)`` returns the force (accel * mass) as the same shape as ``q0``; ``p`` is
    the momentum ``mass * velocity``. Returns ``(qs, ps)``: the ``n_steps + 1`` position
    and momentum states (each a list). Second-order accurate and symplectic, so total
    energy stays bounded over long runs.
    """
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    vec = isinstance(q0, (list, tuple))
    q = _vec(q0)
    p = _vec(p0)
    d = len(q)
    m = mass if isinstance(mass, (list, tuple)) else [mass] * d

    def F(qv):
        return _vec(force(qv if vec else qv[0]))

    qs = [list(q)]
    ps = [list(p)]
    f = F(q)
    for _ in range(n_steps):
        # p half-step, q full-step, recompute force, p half-step.
        p_half = [p[i] + 0.5 * dt * f[i] for i in range(d)]
        q = [q[i] + dt * p_half[i] / m[i] for i in range(d)]
        f = F(q)
        p = [p_half[i] + 0.5 * dt * f[i] for i in range(d)]
        qs.append(list(q))
        ps.append(list(p))
    return qs, ps


def leapfrog(force, q0, v0, dt, n_steps):
    """Leapfrog (kick-drift-kick) integration of ``q'' = force(q)`` (unit mass, velocity form).

    ``force(q)`` returns the acceleration. Returns ``(qs, vs)`` positions and velocities.
    Algebraically equivalent to :func:`velocity_verlet` with unit mass; kept as the
    velocity-space form common in N-body simulation.
    """
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    vec = isinstance(q0, (list, tuple))
    q = _vec(q0)
    v = _vec(v0)
    d = len(q)

    def A(qv):
        return _vec(force(qv if vec else qv[0]))

    qs = [list(q)]
    vs = [list(v)]
    a = A(q)
    for _ in range(n_steps):
        v_half = [v[i] + 0.5 * dt * a[i] for i in range(d)]
        q = [q[i] + dt * v_half[i] for i in range(d)]
        a = A(q)
        v = [v_half[i] + 0.5 * dt * a[i] for i in range(d)]
        qs.append(list(q))
        vs.append(list(v))
    return qs, vs
