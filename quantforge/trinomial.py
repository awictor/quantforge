"""Trinomial lattice for American options, with Richardson extrapolation.

A trinomial tree lets the asset move up, stay, or move down at each step. For
the same number of steps it is smoother and converges faster than the
Cox-Ross-Rubinstein binomial tree, because the extra middle branch reduces the
oscillation of the price in the number of steps.

``trinomial_price`` prices an American (or European) option on the Boyle (1986)
lattice. ``richardson_american`` combines two trinomial solves at ``n`` and
``2n`` steps to cancel the leading O(1/n) error term, giving a materially more
accurate price for the same rough cost.
"""

import math

from .bsm import OptionType, _coerce_type, _validate


def trinomial_price(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                    steps=200, american=True):
    """Price an option on a Boyle trinomial lattice.

    Args:
        american: if True, allow early exercise at every node; if False, price
            the European payoff (useful as a convergence cross-check).
        b: cost of carry (defaults to r). Dividend yield q enters as b = r - q.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if steps < 1:
        raise ValueError("steps must be >= 1")

    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        disc = math.exp(-r * t)
        payoff = max(fwd - K, 0.0) if ot is OptionType.CALL else max(K - fwd, 0.0)
        intrinsic = max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)
        return max(disc * payoff, intrinsic) if american else disc * payoff

    dt = t / steps
    # Boyle: space step lambda*sigma*sqrt(dt) with lambda=sqrt(3) is optimal.
    dx = sigma * math.sqrt(3.0 * dt)
    u = math.exp(dx)
    d = 1.0 / u

    # Risk-neutral up/mid/down probabilities under carry b.
    nu = b - 0.5 * sigma * sigma
    pu = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx)
                + nu * dt / dx)
    pd = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx)
                - nu * dt / dx)
    pm = 1.0 - pu - pd
    if min(pu, pm, pd) < 0.0:
        raise ValueError(
            f"unstable tree probabilities (pu={pu:.3g}, pm={pm:.3g}, pd={pd:.3g}); "
            "increase steps"
        )
    disc = math.exp(-r * dt)

    sign = 1.0 if ot is OptionType.CALL else -1.0

    # Terminal layer has 2*steps+1 nodes, index j=0..2*steps maps to
    # price S * u^(steps - j).
    n_nodes = 2 * steps + 1
    values = [0.0] * n_nodes
    for j in range(n_nodes):
        price_node = S * (u ** (steps - j))
        values[j] = max(sign * (price_node - K), 0.0)

    # Backward induction. At layer m the live nodes are j=0..2*m, whose price is
    # S * u^(m - j).
    for m in range(steps - 1, -1, -1):
        for j in range(2 * m + 1):
            cont = disc * (pu * values[j] + pm * values[j + 1] + pd * values[j + 2])
            if american:
                price_node = S * (u ** (m - j))
                exercise = max(sign * (price_node - K), 0.0)
                values[j] = max(cont, exercise)
            else:
                values[j] = cont
    return values[0]


def richardson_american(S, K, t, r, sigma, option_type=OptionType.CALL, b=None,
                        steps=200):
    """Richardson-extrapolated American price from ``n`` and ``2n`` trinomial solves.

    The trinomial price converges to the true value with a leading O(1/n) error,
    so ``2 * P(2n) - P(n)`` cancels that term and converges faster. Returns the
    extrapolated price.
    """
    p_n = trinomial_price(S, K, t, r, sigma, option_type, b, steps=steps, american=True)
    p_2n = trinomial_price(S, K, t, r, sigma, option_type, b, steps=2 * steps, american=True)
    return 2.0 * p_2n - p_n
