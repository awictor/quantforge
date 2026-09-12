"""Installment options on a binomial tree.

An installment option is paid for in a stream of premiums rather than a single
upfront cost: the buyer pays an initial amount, then a fixed installment ``q`` at
each of several dates before expiry. At every installment date the holder may
*lapse* -- refuse to pay and walk away, forfeiting the option for zero value --
so the option is kept alive only while its continuation value exceeds the next
installment. This optionality to abandon makes it a compound option, priced here
by backward induction on a Cox-Ross-Rubinstein tree.

``installment_call`` returns the *fair upfront* value: the amount a buyer should
pay at inception given the agreed installment schedule. With no installments it
reduces to a plain European call on the same tree. A larger installment lowers the
upfront value (more of the cost is deferred, and the holder lapses more often); a
prohibitively large installment drives it to zero.

Pure standard library.
"""

import math


def installment_call(S, K, t, r, sigma, installment, pay_times, steps=200, q=0.0):
    """Fair upfront value of a European installment call on a CRR tree.

    Parameters
    ----------
    S, K, t, r, sigma : float
        Spot, strike, time to expiry (years), risk-free rate, volatility.
    installment : float
        Amount paid at each date in ``pay_times`` to keep the option alive.
        Must be non-negative.
    pay_times : sequence of float
        Installment payment dates in years, each in ``(0, t)``. Duplicates and
        out-of-range values are rejected.
    steps : int
        Number of tree steps.
    q : float
        Continuous dividend yield (carry ``b = r - q``).

    Returns
    -------
    float
        Fair upfront premium at inception. At least zero; equals the plain
        European call value when ``installment`` is zero or ``pay_times`` empty.
    """
    if S <= 0 or K <= 0 or t <= 0 or sigma <= 0 or steps < 1:
        raise ValueError("S, K, t, sigma must be positive and steps >= 1")
    if installment < 0:
        raise ValueError("installment must be non-negative")
    pay_times = list(pay_times)
    if any(not (0.0 < pt < t) for pt in pay_times):
        raise ValueError("every pay_time must lie strictly inside (0, t)")
    if len(set(pay_times)) != len(pay_times):
        raise ValueError("pay_times must be distinct")

    dt = t / steps
    u = math.exp(sigma * math.sqrt(dt))
    d = 1.0 / u
    disc = math.exp(-r * dt)
    p = (math.exp((r - q) * dt) - d) / (u - d)
    if not (0.0 < p < 1.0):
        raise ValueError("risk-neutral probability out of (0, 1)")

    # Map each installment date to the nearest tree step (>0, <steps). If two
    # dates round to the same step, the installments accrue at that step.
    pay_at_step = {}
    for pt in pay_times:
        i = round(pt / dt)
        i = min(max(i, 1), steps - 1)
        pay_at_step[i] = pay_at_step.get(i, 0.0) + installment

    # Terminal call payoff (no installment due at expiry).
    values = []
    for j in range(steps + 1):
        s = S * (u ** j) * (d ** (steps - j))
        values.append(max(s - K, 0.0))

    for i in range(steps - 1, -1, -1):
        due = pay_at_step.get(i, 0.0)
        new = []
        for j in range(i + 1):
            cont = disc * (p * values[j + 1] + (1.0 - p) * values[j])
            if due > 0.0:
                # Holder pays the installment only if worth it; else lapses to 0.
                cont = max(cont - due, 0.0)
            new.append(cont)
        values = new
    return values[0]
