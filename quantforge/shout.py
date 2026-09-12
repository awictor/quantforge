"""Shout options on a binomial tree.

A shout call lets the holder "shout" once before expiry to lock in the current
intrinsic value ``S_shout - K`` as a guaranteed minimum payoff, while keeping the
upside: the terminal payoff is ``max(S_T - K, S_shout - K)``. It is priced by
backward induction on a Cox-Ross-Rubinstein tree, where at each node the holder
takes the max of continuing (no shout yet) or shouting (locking the intrinsic and
holding a residual call). A shout is worth at least the vanilla. Pure standard
library.
"""

import math

from .bsm import call_price, put_price


def shout_call(S, K, t, r, sigma, steps=200, q=0.0):
    """Shout call price on a CRR tree.

    Once shouted at spot ``S*``, the remaining claim pays
    ``max(S_T - K, S* - K)`` -- a guaranteed ``S* - K`` plus a call struck at ``S*``.
    Its value at the shout node is ``(S* - K) e^{-r tau} + call(S*, S*, tau)`` for
    remaining time ``tau``, taken only when positive. Backward induction compares
    shouting versus continuing. At least the vanilla call value.
    """
    if S <= 0 or K <= 0 or t <= 0 or sigma <= 0 or steps < 1:
        raise ValueError("S, K, t, sigma must be positive and steps >= 1")
    dt = t / steps
    u = math.exp(sigma * math.sqrt(dt))
    d = 1.0 / u
    disc = math.exp(-r * dt)
    p = (math.exp((r - q) * dt) - d) / (u - d)
    if not (0.0 < p < 1.0):
        raise ValueError("risk-neutral probability out of (0, 1)")

    def shout_value(s_node, tau):
        # Value of shouting now at spot s_node with remaining time tau.
        intrinsic = s_node - K
        if intrinsic <= 0.0:
            return 0.0
        # Locked intrinsic (paid at expiry) plus an at-the-money-forward call.
        return intrinsic * math.exp(-r * tau) + call_price(s_node, s_node, tau, r,
                                                           sigma, b=r - q)

    # Terminal payoffs (no shout used) = vanilla call intrinsic.
    values = []
    for j in range(steps + 1):
        s = S * (u ** j) * (d ** (steps - j))
        values.append(max(s - K, 0.0))
    for i in range(steps - 1, -1, -1):
        new = []
        tau = (steps - i) * dt
        for j in range(i + 1):
            s = S * (u ** j) * (d ** (i - j))
            cont = disc * (p * values[j + 1] + (1.0 - p) * values[j])
            new.append(max(cont, shout_value(s, tau)))
        values = new
    return values[0]
