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


def _cn_solve(S, K, t, r, sigma, ot, b, american, local_vol_fn,
              n_space, n_time, s_max_mult, psor_tol, psor_max_iter, dt_override):
    """Run the CN march; return (grid, ds, V, V_prev, dt).

    ``V`` is the value grid at ``t=0``; ``V_prev`` is the grid one time step
    earlier (i.e. at ``dt`` remaining flipped -- used for a theta difference).
    """
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
    V_prev = None   # value grid one step before the last (for theta)

    def vol_at(Si, tau):
        # tau = time remaining; forward time = t - tau.
        if local_vol_fn is not None:
            return local_vol_fn(Si, t - tau)
        return sigma

    # March backward in time (tau = time to maturity increases each step).
    for n in range(n_time):
        if n == n_time - 1:
            V_prev = list(V)   # snapshot at one step (dt) remaining
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

    return grid, ds, V, V_prev, dt


def _interp(grid, ds, V, S):
    n_space = len(grid) - 1
    j = min(int(S / ds), n_space - 1)
    w = (S - grid[j]) / ds
    return (1.0 - w) * V[j] + w * V[j + 1]


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

    grid, ds, V, _Vp, _dt = _cn_solve(S, K, t, r, sigma, ot, b, american,
                                      local_vol_fn, n_space, n_time, s_max_mult,
                                      psor_tol, psor_max_iter, None)
    return _interp(grid, ds, V, S)


def crank_nicolson_greeks(S, K, t, r, sigma=None, option_type=OptionType.CALL,
                          b=None, american=False, local_vol_fn=None,
                          n_space=200, n_time=200, s_max_mult=4.0,
                          psor_tol=1e-8, psor_max_iter=10000):
    """Price plus delta, gamma and theta read straight off the CN grid.

    Delta and gamma come from central finite differences of the final value
    grid in spot (no extra solves), and theta from the difference between the
    ``t=0`` grid and the grid one time step earlier. Returns a dict with price,
    delta, gamma and theta (calendar, per year).
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if b is None:
        b = r
    if sigma is None and local_vol_fn is None:
        raise ValueError("provide sigma or local_vol_fn")
    if t == 0:
        raise ValueError("t must be positive for PDE Greeks")

    grid, ds, V, V_prev, dt = _cn_solve(S, K, t, r, sigma, ot, b, american,
                                        local_vol_fn, n_space, n_time,
                                        s_max_mult, psor_tol, psor_max_iter, None)
    n_space_ = len(grid) - 1
    j = min(max(int(S / ds), 1), n_space_ - 2)
    price = _interp(grid, ds, V, S)

    # Node-level central differences, then interpolate to S (which need not land
    # on a grid node) so delta/gamma are unbiased by the grid offset.
    def node_delta(i):
        return (V[i + 1] - V[i - 1]) / (2.0 * ds)

    def node_gamma(i):
        return (V[i + 1] - 2.0 * V[i] + V[i - 1]) / (ds * ds)

    w = (S - grid[j]) / ds
    delta = (1.0 - w) * node_delta(j) + w * node_delta(j + 1)
    gamma = (1.0 - w) * node_gamma(j) + w * node_gamma(j + 1)
    # Calendar theta = -d(value)/d(time-to-maturity), interpolated to S. V has
    # the full maturity, V_prev one step (dt) less.
    if V_prev is not None:
        th_j = -(V[j] - V_prev[j]) / dt
        th_j1 = -(V[j + 1] - V_prev[j + 1]) / dt
        theta = (1.0 - w) * th_j + w * th_j1
    else:
        theta = 0.0
    return {"price": price, "delta": delta, "gamma": gamma, "theta": theta}


def crank_nicolson_barrier(S, K, H, t, r, sigma=None,
                           option_type=OptionType.CALL, barrier="down-out",
                           b=None, rebate=0.0, local_vol_fn=None,
                           n_space=400, n_time=400):
    """Price a continuously-monitored single-barrier option by a CN PDE.

    Knock-out barriers are imposed as an *absorbing* boundary: at every time
    step the value is set to ``rebate`` on the dead side of ``H`` (V = 0 there
    for a zero rebate), the exact continuous-monitoring condition. ``barrier`` is
    ``down-out``/``up-out`` for the directly-solved knock-outs, or
    ``down-in``/``up-in``, obtained from in + out = vanilla (same rebate handling
    as :func:`quantforge.barrier_option`).

    Supports a constant ``sigma`` or a ``local_vol_fn(S, t)`` and a carry
    ``b = r - q``. Returns the value at spot ``S``.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0 or H <= 0:
        raise ValueError("S, K and H must be positive")
    if b is None:
        b = r
    if sigma is None and local_vol_fn is None:
        raise ValueError("provide sigma or local_vol_fn")
    barrier = str(barrier).lower()
    if barrier not in ("down-out", "up-out", "down-in", "up-in"):
        raise ValueError("barrier must be down-out/up-out/down-in/up-in")

    # Knock-in via in + out = vanilla (both with the same rebate at expiry).
    if barrier.endswith("in"):
        out_kind = barrier.replace("in", "out")
        ko = crank_nicolson_barrier(S, K, H, t, r, sigma, ot, out_kind, b=b,
                                    rebate=0.0, local_vol_fn=local_vol_fn,
                                    n_space=n_space, n_time=n_time)
        vanilla = crank_nicolson_price(S, K, t, r, sigma, ot, b=b,
                                       local_vol_fn=local_vol_fn,
                                       n_space=n_space, n_time=n_time)
        return vanilla - ko

    up = barrier.startswith("up")
    call = ot is OptionType.CALL
    if t == 0:
        alive = (S > H) if up is False else (S < H)
        # Note: down-out alive if S > H; up-out alive if S < H.
        alive = (S > H) if not up else (S < H)
        if not alive:
            return rebate
        return max(S - K, 0.0) if call else max(K - S, 0.0)

    # Grid: choose ds so a node lands exactly on the barrier H (otherwise the
    # absorbing condition is applied at the nearest node, a level offset from H,
    # which biases the price by O(ds)). Pick a target ds, then set the number of
    # sub-intervals below H so H is a grid point, and extend to S_max above.
    ref_vol = sigma if sigma is not None else local_vol_fn(S, t)
    s_max = max(S, K, H) * 4.0 * max(1.0, math.exp(ref_vol * math.sqrt(t)))
    ds_target = s_max / n_space
    # Put a node exactly on the barrier so the absorbing V = rebate condition is
    # applied at H itself (not the nearest node). The barrier discretisation is
    # then O(ds); use a fine grid (n_space >= ~1000) for tight agreement with the
    # continuous-monitoring closed form.
    n_below = max(1, int(round(H / ds_target)))
    ds = H / n_below                     # H = n_below * ds exactly
    n_space = max(int(round(s_max / ds)), n_below + 2)
    grid = [i * ds for i in range(n_space + 1)]
    dt = t / n_time

    def payoff(Si):
        return max(Si - K, 0.0) if call else max(K - Si, 0.0)

    def alive_node(Si):
        return (Si > H) if not up else (Si < H)

    V = [payoff(Si) if alive_node(Si) else rebate for Si in grid]

    def vol_at(Si, tau):
        return local_vol_fn(Si, t - tau) if local_vol_fn is not None else sigma

    for n in range(n_time):
        tau_new = (n + 1) * dt
        sub = [0.0] * (n_space + 1)
        diag = [0.0] * (n_space + 1)
        sup = [0.0] * (n_space + 1)
        rhs = [0.0] * (n_space + 1)
        for i in range(1, n_space):
            Si = grid[i]
            vol = vol_at(Si, tau_new)
            sig2 = vol * vol * Si * Si / (ds * ds)
            drift = b * Si / (2.0 * ds)
            alpha = 0.5 * (0.5 * sig2 - drift)
            gamma = 0.5 * (0.5 * sig2 + drift)
            beta = -0.5 * sig2 - 0.5 * r
            sub[i] = -dt * alpha
            diag[i] = 1.0 - dt * beta
            sup[i] = -dt * gamma
            rhs[i] = (V[i] + dt * alpha * V[i - 1] + dt * beta * V[i]
                      + dt * gamma * V[i + 1])
        disc = math.exp(-r * tau_new)
        carry_fac = math.exp((b - r) * tau_new)
        # Live-side far boundary; dead side gets the rebate.
        if call:
            lowV = 0.0
            highV = grid[n_space] * carry_fac - K * disc
        else:
            lowV = K * disc
            highV = 0.0
        diag[0] = 1.0
        sup[0] = 0.0
        rhs[0] = rebate if not alive_node(grid[0]) else lowV
        sub[n_space] = 0.0
        diag[n_space] = 1.0
        rhs[n_space] = rebate if not alive_node(grid[n_space]) else highV

        V = _thomas(sub, diag, sup, rhs)
        # Absorbing barrier: overwrite the dead side with the rebate.
        for i in range(n_space + 1):
            if not alive_node(grid[i]):
                V[i] = rebate

    return _interp(grid, ds, V, S)
