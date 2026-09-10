"""Crank-Nicolson PDE solver for European and American options.

Solves the Black-Scholes-Merton PDE on a spot grid by the Crank-Nicolson scheme
(second-order accurate in time and space, unconditionally stable). A constant
volatility gives vanilla BSM; a ``local_vol_fn(S, t)`` gives a local-volatility
solver. American exercise is handled by projected successive over-relaxation
(PSOR): each implicit solve is iterated with the early-exercise constraint
``V >= payoff`` applied node-by-node, which converges to the linear-
complementarity solution of the free-boundary problem.

The grid is uniform in spot on ``[0, S_max]`` with ``S_max`` several standard
deviations above the strike; the boundary conditions are the standard vanilla
ones (call: ``V -> S - K e^{-r tau}`` at the top, 0 at 0; put: mirror). Pure
standard library.
"""

import math

from .bsm import OptionType, _coerce_type


def _thomas(a, b, c, d):
    """Solve a tridiagonal system (sub a, diag b, super c, rhs d)."""
    n = len(b)
    cp = [0.0] * n
    dp = [0.0] * n
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / m if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = [0.0] * n
    x[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def crank_nicolson_price(S, K, t, r, sigma=None, option_type=OptionType.CALL,
                         b=None, american=False, local_vol_fn=None,
                         n_space=200, n_time=200, s_max_mult=4.0,
                         psor_tol=1e-8, psor_max_iter=10000):
    """Price a European or American option by a Crank-Nicolson PDE solve.

    Provide either a constant ``sigma`` or a ``local_vol_fn(S, t)`` (time ``t``
    measured forward from today). ``b`` is the cost of carry (defaults to ``r``);
    dividend yield ``q`` enters as ``b = r - q``. American exercise uses PSOR.

    Returns the option value interpolated at spot ``S``.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if b is None:
        b = r
    if sigma is None and local_vol_fn is None:
        raise ValueError("provide sigma or local_vol_fn")
    if t == 0:
        return max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)

    call = ot is OptionType.CALL
    # Spot grid.
    ref_vol = sigma if sigma is not None else local_vol_fn(S, t)
    s_max = max(S, K) * s_max_mult * max(1.0, math.exp(ref_vol * math.sqrt(t)))
    ds = s_max / n_space
    grid = [i * ds for i in range(n_space + 1)]
    dt = t / n_time

    def payoff(Si):
        return max(Si - K, 0.0) if call else max(K - Si, 0.0)

    V = [payoff(Si) for Si in grid]

    def vol_at(Si, tau):
        # tau = time remaining; forward time = t - tau.
        if local_vol_fn is not None:
            return local_vol_fn(Si, t - tau)
        return sigma

    # March backward in time (tau = time to maturity increases each step).
    for n in range(n_time):
        tau_new = (n + 1) * dt
        # Interior CN coefficients (node-dependent for local vol).
        sub = [0.0] * (n_space + 1)
        diag = [0.0] * (n_space + 1)
        sup = [0.0] * (n_space + 1)
        rhs = [0.0] * (n_space + 1)
        for i in range(1, n_space):
            Si = grid[i]
            vol = vol_at(Si, tau_new)
            sig2 = vol * vol * Si * Si / (ds * ds)
            drift = b * Si / (2.0 * ds)
            a_i = 0.5 * (sig2 - 0.0)          # diffusion
            # alpha, beta, gamma of the operator L on interior nodes.
            alpha = 0.5 * (0.5 * sig2 - drift)     # coefficient of V[i-1]
            gamma = 0.5 * (0.5 * sig2 + drift)     # coefficient of V[i+1]
            beta = -0.5 * sig2 - 0.5 * r           # coefficient of V[i]
            # CN: (I - dt/1 * theta L) V_new = (I + dt*(1-theta) L) V_old, theta=0.5
            sub[i] = -dt * alpha
            diag[i] = 1.0 - dt * beta
            sup[i] = -dt * gamma
            rhs[i] = (V[i]
                      + dt * alpha * V[i - 1]
                      + dt * beta * V[i]
                      + dt * gamma * V[i + 1])
        # Boundary conditions.
        disc = math.exp(-r * tau_new)
        carry_fac = math.exp((b - r) * tau_new)
        if call:
            lowV = 0.0
            highV = grid[n_space] * carry_fac - K * disc
        else:
            lowV = K * disc
            highV = 0.0
        diag[0] = 1.0
        sup[0] = 0.0
        rhs[0] = lowV
        sub[n_space] = 0.0
        diag[n_space] = 1.0
        rhs[n_space] = highV

        if not american:
            V = _thomas(sub, diag, sup, rhs)
        else:
            # PSOR with the early-exercise constraint V >= payoff.
            payoff_grid = [payoff(Si) for Si in grid]
            V[0] = max(lowV, payoff_grid[0])
            V[n_space] = max(highV, payoff_grid[n_space])
            omega = 1.5
            for _ in range(psor_max_iter):
                err = 0.0
                for i in range(1, n_space):
                    resid = (rhs[i] - sub[i] * V[i - 1] - sup[i] * V[i + 1]) \
                        / diag[i]
                    new = max(payoff_grid[i], V[i] + omega * (resid - V[i]))
                    err = max(err, abs(new - V[i]))
                    V[i] = new
                if err < psor_tol:
                    break

    # Interpolate the value at S.
    j = min(int(S / ds), n_space - 1)
    w = (S - grid[j]) / ds
    return (1.0 - w) * V[j] + w * V[j + 1]
