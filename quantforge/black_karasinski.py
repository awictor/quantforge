"""Black-Karasinski short-rate model (log-normal, trinomial tree).

Black-Karasinski models the *log* short rate as a mean-reverting
Ornstein-Uhlenbeck process,

    d ln r = kappa (theta - ln r) dt + sigma dW,

so the rate ``r = exp(ln r)`` is strictly positive -- no negative rates, unlike
Vasicek/Hull-White. There is no closed-form bond price, so it is priced on a
Hull-White-style trinomial tree in ``x = ln r``: build a mean-reverting trinomial
lattice for ``x``, exponentiate to get the rate at each node, and discount a unit
cashflow back with backward induction. Pure standard library.
"""

import math


def bk_zero_coupon_bond(r0, kappa, theta, sigma, t, steps=50):
    """Zero-coupon bond price under Black-Karasinski via a trinomial tree.

    Parameters
    ----------
    r0 : float
        Current short rate (> 0).
    kappa : float
        Mean-reversion speed of the log rate (> 0).
    theta : float
        Long-run mean of the *log* rate (``ln`` of the target rate level).
    sigma : float
        Volatility of the log rate (> 0).
    t : float
        Bond maturity (years).
    steps : int
        Number of tree steps.

    Returns
    -------
    float
        Price of a unit zero-coupon bond maturing at ``t``. Strictly in ``(0, 1]``
        for positive rates; falls as ``r0`` or ``sigma`` rises.
    """
    if r0 <= 0 or kappa <= 0 or sigma <= 0 or t <= 0 or steps < 1:
        raise ValueError("r0, kappa, sigma, t must be positive and steps >= 1")

    dt = t / steps
    dx = sigma * math.sqrt(3.0 * dt)        # Hull-White standard spacing
    x0 = math.log(r0)
    m = -kappa * dt                         # mean-reversion drift factor on (x - theta)

    # Trinomial branching probabilities (Hull-White), constant here since the
    # drift is linear and spacing uniform; computed at each node from its level.
    def probs(j):
        # Expected change in x over dt from node level x = x0 + j*dx (approx via
        # deviation from theta). Use standard HW trinomial with mean reversion.
        alpha = m * (x0 + j * dx - theta)   # E[dx] contribution
        eta = alpha / dx
        pu = 1.0 / 6.0 + 0.5 * (eta * eta + eta)
        pm = 2.0 / 3.0 - eta * eta
        pd = 1.0 / 6.0 + 0.5 * (eta * eta - eta)
        return pu, pm, pd

    # Node index range grows by 1 each step; cap width for stability.
    # Value at maturity = 1 for every node.
    # Represent values by absolute index j in [-steps, steps].
    width = steps
    values = {j: 1.0 for j in range(-width, width + 1)}

    for _ in range(steps):
        new = {}
        for j in range(-width, width + 1):
            x = x0 + j * dx
            r = math.exp(x)
            disc = math.exp(-r * dt)
            pu, pm, pd = probs(j)
            up = values.get(j + 1, values[max(-width, min(width, j + 1))])
            mid = values[j]
            dn = values.get(j - 1, values[max(-width, min(width, j - 1))])
            new[j] = disc * (pu * up + pm * mid + pd * dn)
        values = new
    return values[0]
