"""Option-position VaR/ES under the smile-implied risk-neutral density.

Given a position whose horizon P&L is a function of the terminal spot,
``pnl(S_T)``, and the option smile, the risk-neutral density
(:mod:`quantforge.rnd`) gives the full P&L distribution -- so its Value-at-Risk
and Expected Shortfall follow by integrating the density.

This is a *risk-neutral* (Q-measure) VaR: useful for a model-free, smile-
consistent tail read on a structured payoff, though it is not the physical-
measure regulatory VaR. VaR is returned as a positive loss number; ES is the
mean loss conditional on breaching it. Pure standard library.
"""

import math

from .rnd import density_grid_from_smile


def density_var_es(S0, t, r, vol_fn, pnl, confidence=0.99, q=0.0, n=800,
                   width=8.0):
    """Risk-neutral VaR and Expected Shortfall of a payoff ``pnl(S_T)``.

    Args:
        pnl: horizon profit-and-loss as a function of the terminal spot (losses
            negative).
        confidence: e.g. 0.99 for the 99% level.

    Returns ``(var, es)`` with both as positive loss numbers: ``var`` is the loss
    the P&L does not exceed with probability ``confidence`` under the
    risk-neutral density, and ``es`` the mean loss beyond it.
    """
    ks, dens = density_grid_from_smile(S0, t, r, vol_fn, q=q, n=n, width=width)
    dens = [max(d, 0.0) for d in dens]
    dK = ks[1] - ks[0]
    mass = sum(dens) * dK
    if mass <= 0:
        raise ValueError("degenerate density")
    dens = [d / mass for d in dens]

    # Build (pnl, prob) atoms and sort by P&L ascending (worst first).
    atoms = [(pnl(ks[i]), dens[i] * dK) for i in range(len(ks))]
    atoms.sort(key=lambda a: a[0])

    alpha = 1.0 - confidence          # tail mass
    cum = 0.0
    var_pnl = atoms[0][0]
    for pl, p in atoms:
        cum += p
        if cum >= alpha:
            var_pnl = pl              # the alpha-quantile of P&L
            break

    # ES: mean P&L over the worst alpha mass (up to and including the quantile).
    cum = 0.0
    num = 0.0
    for pl, p in atoms:
        remaining = alpha - cum
        if remaining <= 0:
            break               # filled the tail mass
        take = min(p, remaining)
        num += pl * take        # p may be ~0 in the deep tail; that is fine
        cum += take
    es_pnl = num / alpha if alpha > 0 else var_pnl

    return -var_pnl, -es_pnl          # report as positive losses
