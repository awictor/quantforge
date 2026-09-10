"""Leisen-Reimer binomial tree for European and American options.

The Cox-Ross-Rubinstein tree converges slowly and its error oscillates with the
step count (the strike drifts between adjacent tree nodes). Leisen & Reimer
(1996) instead centre the tree on the strike and set the up-move and
probability from a Peizer-Pratt normal-approximation *inversion* of the
Black-Scholes ``d1``/``d2``. The result converges at ``O(1/n^2)`` and, crucially,
monotonically -- so a modest odd step count (a few dozen) matches what CRR needs
thousands of steps to reach, and Richardson extrapolation on two sizes is clean.

For European payoffs the tree reproduces Black-Scholes; for American exercise it
adds the usual early-exercise check at every node. Pure standard library.
"""

import math

from .bsm import OptionType, _coerce_type, _validate


def _peizer_pratt(z, n):
    """Peizer-Pratt inversion: a normal-CDF approximation used as a tree prob.

    Returns ``h(z)`` in (0, 1), the discrete probability that reproduces
    ``N(z)`` on an ``n``-step tree (method 2 of Leisen-Reimer). ``n`` must be
    odd so the strike sits at the centre of the terminal nodes.
    """
    c = 1.0 if z >= 0 else -1.0
    inner = z / (n + 1.0 / 3.0 + 0.1 / (n + 1.0))
    return 0.5 + c * 0.5 * math.sqrt(1.0 - math.exp(-inner * inner * (n + 1.0 / 6.0)))


def leisen_reimer_price(S, K, t, r, sigma, option_type=OptionType.CALL,
                        b=None, steps=101, american=False) -> float:
    """Price an option with a Leisen-Reimer binomial tree.

    Args:
        steps: number of time steps; forced to the next odd integer so the tree
            straddles the strike (the source of the fast, monotone convergence).
        american: if ``True``, apply an early-exercise check at each node;
            otherwise price the European payoff.
        b: cost of carry (defaults to ``r``); dividend yield ``q`` enters as
            ``b = r - q``.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if steps < 1:
        raise ValueError("steps must be >= 1")
    n = steps if steps % 2 == 1 else steps + 1  # odd step count

    if t == 0 or sigma == 0:
        fwd = S * math.exp(b * t)
        disc = math.exp(-r * t)
        payoff = max(fwd - K, 0.0) if ot is OptionType.CALL else max(K - fwd, 0.0)
        intrinsic = max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)
        return max(disc * payoff, intrinsic) if american else disc * payoff

    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log(S / K) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    d2 = d1 - vsqrt

    p = _peizer_pratt(d2, n)          # risk-neutral probability
    p_star = _peizer_pratt(d1, n)     # probability under the share measure
    # Up/down moves chosen so the tree is centred on the strike.
    growth = math.exp(b * t / n)
    u = growth * p_star / p
    dd = growth * (1.0 - p_star) / (1.0 - p)
    disc = math.exp(-r * t / n)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    # Terminal payoffs.
    values = []
    for i in range(n + 1):
        price_node = S * (u ** (n - i)) * (dd ** i)
        values.append(max(sign * (price_node - K), 0.0))

    for step in range(n - 1, -1, -1):
        for i in range(step + 1):
            cont = disc * (p * values[i] + (1.0 - p) * values[i + 1])
            if american:
                price_node = S * (u ** (step - i)) * (dd ** i)
                exercise = max(sign * (price_node - K), 0.0)
                values[i] = max(cont, exercise)
            else:
                values[i] = cont

    return values[0]


def leisen_reimer_american_accel(S, K, t, r, sigma, option_type=OptionType.CALL,
                                 b=None, steps=101) -> float:
    """Richardson-extrapolated Leisen-Reimer American price (Broadie-Detemple).

    American LR convergence is only ``O(1/n)`` (the smooth-payoff assumption
    behind the Peizer-Pratt inversion breaks at the early-exercise boundary),
    unlike the ``O(1/n^2)`` European case. Broadie & Detemple (1996) cancel that
    leading ``1/n`` term with a two-point Richardson extrapolation between an
    ``n``-step and a ``2n``-step tree:

        V_ext = 2 * V(2n) - V(n).

    For the same work this is several times more accurate than a single tree, so
    a moderate ``steps`` reaches four-figure accuracy. ``b`` is the cost of carry
    (dividend yield ``q`` via ``b = r - q``).
    """
    vn = leisen_reimer_price(S, K, t, r, sigma, option_type, b=b,
                             steps=steps, american=True)
    v2n = leisen_reimer_price(S, K, t, r, sigma, option_type, b=b,
                              steps=2 * steps, american=True)
    return 2.0 * v2n - vn


def leisen_reimer_greeks(S, K, t, r, sigma, option_type=OptionType.CALL,
                         b=None, steps=101, american=False):
    """Delta and gamma from Leisen-Reimer tree nodes, plus FD vega/theta.

    Delta and gamma are read directly off the first two time steps of the tree
    (no extra pricing passes), while vega and theta use small central
    differences. Returns a dict with price, delta, gamma, vega and theta.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        return leisen_reimer_price(S_, K, t_, r, sigma_, ot, b=b,
                                   steps=steps, american=american)

    base = px()
    hS = 1e-3 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.5 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma,
            "vega": vega, "theta": theta}
