"""Arithmetic-average Asian option by an augmented-state PDE.

A continuously-monitored arithmetic Asian depends on the running integral
``I_t = int_0^t S_u du`` as well as the spot, so its value solves a 2D PDE in
``(S, I)``:

    V_t + 0.5 sigma^2 S^2 V_SS + b S V_S + S V_I - r V = 0,

with terminal payoff ``(I_T / T - K)^+`` for a fixed-strike average-price call.
The integral state does not diffuse -- ``dI = S dt`` is pure transport -- so the
step is operator-split: a semi-Lagrangian transport in ``I`` (each backward
step, ``I`` moves by ``S dt`` along its characteristic, handled by
interpolation) followed by a Crank-Nicolson diffusion in ``S`` at each fixed
``I`` line. Cross-checked against the Turnbull-Wakeman closed form and the
arithmetic-Asian Monte Carlo. Pure standard library.
"""

import math

from .pde import _thomas
from .bsm import OptionType, _coerce_type


def asian_pde_price(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                    n_s=120, n_i=120, n_time=100, s_max_mult=4.0):
    """Fixed-strike continuously-averaged arithmetic Asian by a 2D PDE.

    ``b`` is the cost of carry (defaults to ``r``); dividend yield ``q`` enters
    as ``b = r - q``. Averaging runs over the full life ``[0, t]``; the payoff is
    ``(A_t - K)^+`` for a call and ``(K - A_t)^+`` for a put, with
    ``A_t = (1/t) int_0^t S_u du``.

    Returns the value at spot ``S`` (the running integral starts at 0).
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if b is None:
        b = r
    call = ot is OptionType.CALL
    if t == 0:
        return max(S - K, 0.0) if call else max(K - S, 0.0)

    s_max = S * s_max_mult * max(1.0, math.exp(sigma * math.sqrt(t)))
    # The average A = I/t rarely exceeds a few times max(S, K); bounding I at
    # a_max * t (rather than s_max * t) concentrates the I-grid where the payoff
    # kink lives, so the same n_i resolves the strike far better.
    a_max = max(S, K) * (1.0 + 1.5 * max(0.3, sigma) * math.sqrt(t)) * 1.6
    i_max = a_max * t
    ds = s_max / n_s
    di = i_max / n_i
    dt = t / n_time
    s_grid = [k * ds for k in range(n_s + 1)]
    i_grid = [k * di for k in range(n_i + 1)]

    def payoff(Si, Ii):
        A = Ii / t
        return max(A - K, 0.0) if call else max(K - A, 0.0)

    # V[si][ii] terminal.
    V = [[payoff(s_grid[si], i_grid[ii]) for ii in range(n_i + 1)]
         for si in range(n_s + 1)]

    # CN coefficients in S (constant per S node, applied on each I-line).
    disc = math.exp(-r * dt)
    for _ in range(n_time):
        # --- Transport in I: backward step moves I by -S dt along dI = S dt. ---
        # V_transport[si][ii] = V(., I + S*dt) by linear interpolation in I.
        Vt = [[0.0] * (n_i + 1) for _ in range(n_s + 1)]
        for si in range(n_s + 1):
            shift = s_grid[si] * dt
            for ii in range(n_i + 1):
                I_target = i_grid[ii] + shift
                pos = I_target / di
                lo = int(pos)
                if lo >= n_i:
                    Vt[si][ii] = V[si][n_i]
                else:
                    w = pos - lo
                    Vt[si][ii] = (1.0 - w) * V[si][lo] + w * V[si][lo + 1]

        # --- CN diffusion in S on each I-line. ---
        Vn = [[0.0] * (n_i + 1) for _ in range(n_s + 1)]
        for ii in range(n_i + 1):
            sub = [0.0] * (n_s + 1)
            dia = [0.0] * (n_s + 1)
            sup = [0.0] * (n_s + 1)
            rhs = [0.0] * (n_s + 1)
            # S = 0 boundary: pure discount of the transported value.
            dia[0] = 1.0
            rhs[0] = Vt[0][ii] * disc
            # S = s_max boundary: linear in S (gamma -> 0).
            dia[n_s] = 1.0
            rhs[n_s] = Vt[n_s][ii]
            for si in range(1, n_s):
                Ssi = s_grid[si]
                sig2 = sigma * sigma * Ssi * Ssi / (ds * ds)
                drift = b * Ssi / (2.0 * ds)
                alpha = 0.5 * sig2 - drift
                gamma = 0.5 * sig2 + drift
                beta = -sig2 - r
                th = 0.5
                sub[si] = -th * dt * alpha
                dia[si] = 1.0 - th * dt * beta
                sup[si] = -th * dt * gamma
                rhs[si] = (Vt[si][ii]
                           + (1 - th) * dt * alpha * Vt[si - 1][ii]
                           + (1 - th) * dt * beta * Vt[si][ii]
                           + (1 - th) * dt * gamma * Vt[si + 1][ii])
            col = _thomas(sub, dia, sup, rhs)
            for si in range(n_s + 1):
                Vn[si][ii] = col[si]
        V = Vn

    # Interpolate at (S, I=0).
    si = min(max(int(S / ds), 0), n_s - 1)
    w = (S - s_grid[si]) / ds
    return (1.0 - w) * V[si][0] + w * V[si + 1][0]
