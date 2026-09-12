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


def callable_bond_price_with_spread(face, coupon_rate, maturity, r0, sigma,
                                    spread, freq=1, call_price=None,
                                    put_price=None, p=0.5):
    """Callable-bond price with a constant spread added to every tree rate.

    Shifts the whole short-rate lattice up by ``spread`` before backward
    induction, so a positive spread discounts harder and lowers the price. The
    building block for the option-adjusted spread solve.
    """
    if face <= 0 or maturity <= 0 or freq < 1:
        raise ValueError("face, maturity must be positive and freq >= 1")
    n = int(round(maturity * freq))
    dt = 1.0 / freq
    coupon = face * coupon_rate / freq
    rates = _rate_tree(r0, sigma, n, dt)
    values = [face + coupon for _ in range(n + 1)]
    for i in range(n - 1, -1, -1):
        new = []
        for j in range(i + 1):
            disc = math.exp(-(rates[i][j] + spread) * dt)
            cont = disc * (p * values[j + 1] + (1.0 - p) * values[j]) + coupon
            v = cont
            if call_price is not None:
                v = min(v, call_price + coupon)
            if put_price is not None:
                v = max(v, put_price + coupon)
            new.append(v)
        values = new
    return values[0]


def option_adjusted_spread(market_price, face, coupon_rate, maturity, r0, sigma,
                           freq=1, call_price=None, put_price=None, p=0.5,
                           tol=1e-8, max_iter=100):
    """Option-adjusted spread: constant rate spread repricing the bond to market.

    Bisection on the :func:`callable_bond_price_with_spread` (monotone decreasing
    in the spread) to hit ``market_price``. Positive when the market price is below
    the zero-spread model price. Strips out the embedded option so the spread
    reflects credit/liquidity risk on a like-for-like basis.
    """
    if market_price <= 0:
        raise ValueError("market_price must be positive")

    def px(sp):
        return callable_bond_price_with_spread(face, coupon_rate, maturity, r0,
                                               sigma, sp, freq, call_price,
                                               put_price, p)

    lo, hi = -0.5, 1.0
    p_lo, p_hi = px(lo), px(hi)
    if not (min(p_lo, p_hi) - 1e-6 <= market_price <= max(p_lo, p_hi) + 1e-6):
        raise ValueError("market_price outside the achievable spread range")
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        pm = px(mid)
        if abs(pm - market_price) < tol:
            return mid
        if pm > market_price:   # price too high -> raise spread
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


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
