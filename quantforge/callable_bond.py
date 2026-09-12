"""Callable and puttable bond pricing on a short-rate binomial tree.

Prices a bond with embedded call/put options by backward induction on a
recombining binomial short-rate tree (a simple Black-Derman-Toy-style lattice
with constant volatility). At each node the continuation value is the
discounted risk-neutral average of the two successor values plus the coupon; a
callable bond caps the holder's value at the call price (the issuer redeems when
cheaper), a puttable bond floors it at the put price (the holder redeems when the
bond is worth less). Pure standard library.
"""

import math


def _rate_tree(r0, sigma, n, dt):
    """Recombining short-rate tree ``r(i,j) = r0 exp(sigma (2j - i) sqrt(dt))``.

    Log-normal rates (BDT-style) with constant volatility; node ``j`` at step ``i``
    has ``j`` up-moves.
    """
    tree = []
    for i in range(n + 1):
        row = [r0 * math.exp(sigma * (2 * j - i) * math.sqrt(dt)) for j in range(i + 1)]
        tree.append(row)
    return tree


def callable_bond_price(face, coupon_rate, maturity, r0, sigma, freq=1,
                        call_price=None, put_price=None, p=0.5):
    """Price a (possibly callable/puttable) bond on a binomial short-rate tree.

    ``freq`` coupons per year over ``maturity`` years; the lattice has one step per
    coupon. ``call_price`` caps the value at each pre-maturity node (issuer's call);
    ``put_price`` floors it (holder's put); either may be ``None``. With both
    ``None`` this is the straight bond on the tree. Risk-neutral up-probability
    ``p``. Returns the time-zero price.
    """
    if face <= 0 or maturity <= 0 or freq < 1:
        raise ValueError("face, maturity must be positive and freq >= 1")
    n = int(round(maturity * freq))
    dt = 1.0 / freq
    coupon = face * coupon_rate / freq
    rates = _rate_tree(r0, sigma, n, dt)
    # Terminal: redeem face plus final coupon.
    values = [face + coupon for _ in range(n + 1)]
    for i in range(n - 1, -1, -1):
        new = []
        for j in range(i + 1):
            disc = math.exp(-rates[i][j] * dt)
            cont = disc * (p * values[j + 1] + (1.0 - p) * values[j]) + coupon
            v = cont
            if call_price is not None:
                v = min(v, call_price + coupon)   # issuer calls if cheaper
            if put_price is not None:
                v = max(v, put_price + coupon)    # holder puts if worth less
            new.append(v)
        values = new
    return values[0]


def straight_bond_tree_price(face, coupon_rate, maturity, r0, sigma, freq=1,
                             p=0.5):
    """Straight (option-free) bond on the same tree -- the no-optionality baseline."""
    return callable_bond_price(face, coupon_rate, maturity, r0, sigma, freq,
                               call_price=None, put_price=None, p=p)


def call_option_value(face, coupon_rate, maturity, r0, sigma, call_price,
                      freq=1, p=0.5):
    """Value of the embedded call to the issuer: straight minus callable price.

    ``straight_bond - callable_bond`` -- the value the call option strips from the
    bondholder (non-negative). The callable bond trades cheaper by this amount.
    """
    straight = straight_bond_tree_price(face, coupon_rate, maturity, r0, sigma,
                                        freq, p)
    callable_ = callable_bond_price(face, coupon_rate, maturity, r0, sigma, freq,
                                    call_price=call_price, p=p)
    return straight - callable_
