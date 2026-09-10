"""Cox-Ross-Rubinstein binomial lattice for American-exercise options.

Prices early-exercise options that BSM cannot handle in closed form. The tree
converges to the BSM value for European payoffs as ``steps`` grows, which the
test suite uses as a cross-check.
"""

import math

from .bsm import OptionType, _coerce_type, _validate


def american_price(S, K, t, r, sigma, option_type=OptionType.CALL, b=None, steps=500):
    """Price an American option via a CRR binomial tree.

    Args:
        steps: number of time steps. Higher = more accurate, O(steps^2) work.
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
        return disc * payoff

    dt = t / steps
    u = math.exp(sigma * math.sqrt(dt))
    d = 1.0 / u
    # Risk-neutral up-probability under carry b.
    p = (math.exp(b * dt) - d) / (u - d)
    if not (0.0 <= p <= 1.0):
        raise ValueError(
            "unstable tree: increase steps or check inputs "
            f"(p={p:.4g}); dt too large for given sigma/b"
        )
    disc = math.exp(-r * dt)

    # Terminal asset prices and payoffs.
    sign = 1.0 if ot is OptionType.CALL else -1.0
    values = []
    for i in range(steps + 1):
        price_node = S * (u ** (steps - i)) * (d ** i)
        values.append(max(sign * (price_node - K), 0.0))

    # Backward induction with early-exercise check at each node.
    for step in range(steps - 1, -1, -1):
        for i in range(step + 1):
            cont = disc * (p * values[i] + (1.0 - p) * values[i + 1])
            price_node = S * (u ** (step - i)) * (d ** i)
            exercise = max(sign * (price_node - K), 0.0)
            values[i] = max(cont, exercise)

    return values[0]
