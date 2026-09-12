"""Convertible bond pricing on an equity binomial lattice.

Prices a convertible bond by backward induction on a Cox-Ross-Rubinstein tree of
the issuer's stock. At each node the holder takes the max of continuing to hold
(the discounted risk-neutral continuation plus coupon), converting into
``conversion_ratio`` shares, or being put; the issuer caps the value at the call
price. This captures the American-style optionality that the closed-form
component decomposition (``equity_comp.convertible_bond_value``) misses, while
converging to it for a European, non-callable convertible. Pure standard library.
"""

import math


def convertible_bond_lattice(S, sigma, face, conversion_ratio, coupon_rate,
                             maturity, r, steps=200, call_price=None,
                             put_price=None, credit_spread=0.0, q=0.0):
    """Convertible bond value on a CRR equity tree by backward induction.

    ``S`` current stock, ``conversion_ratio`` shares per bond, coupons at rate
    ``coupon_rate`` on ``face`` spread evenly across the steps. At each node the
    holder value is ``max(continuation, conversion_ratio * S_node, put_price)`` and
    the issuer caps it at ``call_price`` (both optional). Continuation discounts at
    ``r + credit_spread`` (risky). ``q`` is the dividend yield. Returns the
    time-zero convertible price.
    """
    if S <= 0 or sigma <= 0 or maturity <= 0 or steps < 1:
        raise ValueError("S, sigma, maturity must be positive and steps >= 1")
    if conversion_ratio < 0 or face <= 0:
        raise ValueError("conversion_ratio non-negative and face positive")
    dt = maturity / steps
    u = math.exp(sigma * math.sqrt(dt))
    d = 1.0 / u
    disc = math.exp(-(r + credit_spread) * dt)
    p = (math.exp((r - q) * dt) - d) / (u - d)
    if not (0.0 < p < 1.0):
        raise ValueError("risk-neutral probability out of (0,1); check inputs")
    coupon = face * coupon_rate / steps
    # Terminal stock prices and payoffs.
    values = []
    for j in range(steps + 1):
        s_node = S * (u ** j) * (d ** (steps - j))
        redeem = max(face + coupon, conversion_ratio * s_node)
        values.append(redeem)
    for i in range(steps - 1, -1, -1):
        new = []
        for j in range(i + 1):
            s_node = S * (u ** j) * (d ** (i - j))
            cont = disc * (p * values[j + 1] + (1.0 - p) * values[j]) + coupon
            conv = conversion_ratio * s_node
            v = max(cont, conv)
            if put_price is not None:
                v = max(v, put_price + coupon)
            if call_price is not None:
                # Issuer calls, but the holder converts if that beats the call.
                v = min(v, max(call_price + coupon, conv))
            new.append(v)
        values = new
    return values[0]
