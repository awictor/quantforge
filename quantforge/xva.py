"""Counterparty valuation adjustments (CVA / DVA / BCVA).

Unilateral credit valuation adjustment discounts the expected loss from a
counterparty defaulting while the trade is in the money to us:

    CVA = LGD * sum_i EE(t_i) * DF(t_i) * [Q(t_{i-1}) - Q(t_i)]

where ``EE`` is the (positive) expected exposure profile, ``DF`` the discount
factors, ``Q`` the counterparty survival probability, and ``LGD = 1 - recovery``.
DVA is the mirror term on our own default against the negative exposure (a
benefit), and the bilateral adjustment is ``BCVA = CVA - DVA``. The default
buckets ``Q(t_{i-1}) - Q(t_i)`` are the marginal default probabilities over each
grid step. Pure standard library; pairs with :class:`quantforge.SurvivalCurve`.
"""

import math


def _disc_fn(r):
    return r if callable(r) else (lambda t: math.exp(-r * t))


def marginal_default_probs(curve, grid_times):
    """Marginal default probability in each grid bucket ``Q(t_{i-1}) - Q(t_i)``.

    ``grid_times`` are the bucket end points (strictly increasing, positive); the
    first bucket runs from 0. Returns one probability per bucket, each in
    ``[0, 1]`` and summing to ``1 - Q(t_last)`` (the total default probability by
    the horizon).
    """
    if not grid_times:
        raise ValueError("need at least one grid time")
    prev_t = 0.0
    prev_q = 1.0
    out = []
    for t in grid_times:
        if t <= prev_t:
            raise ValueError("grid_times must be strictly increasing and positive")
        q = curve.survival(t)
        out.append(prev_q - q)
        prev_t, prev_q = t, q
    return out


def cva(curve, grid_times, expected_exposure, r, recovery=0.4):
    """Unilateral CVA from an expected-exposure profile and a survival curve.

    ``expected_exposure[i]`` is the positive expected exposure at ``grid_times[i]``
    (the representative exposure over bucket ``i``), discounted by ``DF`` (flat
    rate ``r`` or a callable ``r(t)``) and weighted by the counterparty's marginal
    default probability :func:`marginal_default_probs`. Scaled by
    ``LGD = 1 - recovery``. Non-negative, increasing in exposure and in hazard.
    """
    if len(expected_exposure) != len(grid_times):
        raise ValueError("expected_exposure and grid_times must align")
    if not (0.0 <= recovery <= 1.0):
        raise ValueError("recovery must be in [0, 1]")
    df = _disc_fn(r)
    lgd = 1.0 - recovery
    dq = marginal_default_probs(curve, grid_times)
    total = 0.0
    for t, ee, p in zip(grid_times, expected_exposure, dq):
        if ee < 0.0:
            raise ValueError("expected exposure must be non-negative")
        total += ee * df(t) * p
    return lgd * total


def dva(own_curve, grid_times, negative_expected_exposure, r, recovery=0.4):
    """Debit valuation adjustment: the mirror of :func:`cva` on our own default.

    ``negative_expected_exposure[i]`` is the expected exposure of the counterparty
    to us (our negative exposure, entered as a non-negative magnitude). Weighted by
    *our* marginal default probability from ``own_curve`` and ``LGD``. A benefit to
    us, so it is subtracted from CVA in the bilateral adjustment.
    """
    return cva(own_curve, grid_times, negative_expected_exposure, r, recovery)


def bcva(cpty_curve, own_curve, grid_times, epe, ene, r,
         cpty_recovery=0.4, own_recovery=0.4):
    """Bilateral CVA ``BCVA = CVA - DVA``.

    ``epe`` is the expected positive exposure profile (counterparty default risk)
    and ``ene`` the expected negative exposure profile (our default benefit).
    Returns the net adjustment to the risk-free value; positive when counterparty
    risk dominates.
    """
    return (cva(cpty_curve, grid_times, epe, r, cpty_recovery)
            - dva(own_curve, grid_times, ene, r, own_recovery))
