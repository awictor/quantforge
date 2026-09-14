"""One-dimensional SDE integrators: Euler-Maruyama and Milstein.

Numerically integrate ``dX = a(X, t) dt + b(X, t) dW`` from user-supplied drift ``a`` and
diffusion ``b``. Euler-Maruyama is the basic scheme (strong order 0.5); Milstein adds the
``b b' `` correction term for strong order 1.0, converging faster for state-dependent
diffusions like GBM and CIR. Uses a deterministic LCG + Box-Muller normal stream so paths
are reproducible per seed. Pure standard library.
"""

import math


def _normals(seed, n):
    """Deterministic standard-normal stream (LCG + Box-Muller)."""
    state = seed & 0x7FFFFFFF or 1

    def unif():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return (state + 0.5) / 0x80000000

    out = []
    while len(out) < n:
        u1 = unif()
        u2 = unif()
        r = math.sqrt(-2.0 * math.log(u1))
        out.append(r * math.cos(2.0 * math.pi * u2))
        out.append(r * math.sin(2.0 * math.pi * u2))
    return out[:n]


def euler_maruyama(drift, diffusion, x0, t, n_steps, seed=1234567):
    """Euler-Maruyama path of ``dX = drift(X, t) dt + diffusion(X, t) dW``.

    ``drift`` and ``diffusion`` are callables ``(x, t) -> float``. Returns the path as a
    list of ``n_steps + 1`` states at times ``0, dt, ..., t``. Strong order 0.5.
    """
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    dt = t / n_steps
    sdt = math.sqrt(dt)
    z = _normals(seed, n_steps)
    x = x0
    path = [x0]
    tc = 0.0
    for k in range(n_steps):
        dw = sdt * z[k]
        x = x + drift(x, tc) * dt + diffusion(x, tc) * dw
        tc += dt
        path.append(x)
    return path


def milstein(drift, diffusion, diffusion_prime, x0, t, n_steps, seed=1234567):
    """Milstein path with the ``0.5 b b' (dW^2 - dt)`` correction (strong order 1.0).

    Adds the derivative of the diffusion ``diffusion_prime(x, t) = d b/d x`` to the
    Euler-Maruyama step, giving strong order 1.0 for state-dependent diffusions. Same
    signature otherwise; returns the ``n_steps + 1`` state path.
    """
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    dt = t / n_steps
    sdt = math.sqrt(dt)
    z = _normals(seed, n_steps)
    x = x0
    path = [x0]
    tc = 0.0
    for k in range(n_steps):
        dw = sdt * z[k]
        b = diffusion(x, tc)
        x = (x + drift(x, tc) * dt + b * dw
             + 0.5 * b * diffusion_prime(x, tc) * (dw * dw - dt))
        tc += dt
        path.append(x)
    return path


def gbm_paths(mu, sigma, x0, t, n_steps, n_paths, seed=1234567, scheme="milstein"):
    """Simulate geometric Brownian motion paths (``dS = mu S dt + sigma S dW``).

    Convenience wrapper over :func:`milstein` (default) or :func:`euler_maruyama` with the
    GBM drift and diffusion. Returns a list of ``n_paths`` state paths. The sample mean of
    the terminal value approaches the analytic ``x0 exp(mu t)`` as ``n_paths`` grows.
    """
    drift = lambda x, tc: mu * x
    diffusion = lambda x, tc: sigma * x
    dprime = lambda x, tc: sigma
    paths = []
    for p in range(n_paths):
        s = seed + p * 100003
        if scheme == "milstein":
            paths.append(milstein(drift, diffusion, dprime, x0, t, n_steps, seed=s))
        elif scheme == "euler":
            paths.append(euler_maruyama(drift, diffusion, x0, t, n_steps, seed=s))
        else:
            raise ValueError("scheme must be 'milstein' or 'euler'")
    return paths
