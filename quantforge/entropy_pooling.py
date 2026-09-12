"""Entropy pooling: reweight scenarios to honor a view (Meucci 2008).

Given ``n`` scenarios for a quantity ``x`` with prior probabilities ``q`` (usually
uniform ``1/n``), entropy pooling finds the posterior probabilities ``p`` that
satisfy a view -- here a target mean ``E_p[x] = target`` -- while staying as close
as possible to the prior in the Kullback-Leibler sense:

    minimize  sum_i p_i log(p_i / q_i)   s.t.   sum_i p_i = 1,  sum_i p_i x_i = target.

The Lagrangian solution is an exponential tilt of the prior,

    p_i = q_i exp(lambda x_i) / Z(lambda),

with the single multiplier ``lambda`` chosen so the posterior mean equals the
target. Because ``E_p[x]`` is monotone increasing in ``lambda``, a bisection on
``lambda`` converges to the unique solution. This is the minimal-information way to
impose a view on a scenario set, and the reweighted scenarios feed straight into
any downstream risk or allocation calculation. Pure standard library.
"""

import math


def _weighted_mean(x, p):
    return sum(x[i] * p[i] for i in range(len(x)))


def entropy_pooling_mean(x, target, prior=None, tol=1e-12, max_iter=200):
    """Posterior scenario probabilities matching a target mean, min relative entropy.

    Parameters
    ----------
    x : sequence of float
        Scenario values of the quantity being viewed.
    target : float
        Desired posterior mean ``E_p[x]``. Must lie strictly inside
        ``(min(x), max(x))`` -- a mean outside the scenario range is infeasible.
    prior : sequence of float, optional
        Prior probabilities (non-negative, summing to 1). Defaults to uniform.
    tol, max_iter : float, int
        Bisection tolerance on the achieved mean and its iteration cap.

    Returns
    -------
    list[float]
        Posterior probabilities ``p`` (an exponential tilt of the prior). They sum
        to 1, reproduce ``target`` as their mean, and reduce to the prior when
        ``target`` equals the prior mean.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least 2 scenarios")
    if prior is None:
        q = [1.0 / n] * n
    else:
        if len(prior) != n:
            raise ValueError("prior must match x in length")
        s = sum(prior)
        if s <= 0 or any(pi < 0 for pi in prior):
            raise ValueError("prior must be non-negative and sum to a positive value")
        q = [pi / s for pi in prior]

    lo_x, hi_x = min(x), max(x)
    if not (lo_x < target < hi_x):
        raise ValueError("target mean must lie strictly inside (min(x), max(x))")

    def posterior(lam):
        # Numerically stable softmax-style tilt.
        mx = max(lam * xi for xi in x)
        w = [q[i] * math.exp(lam * x[i] - mx) for i in range(n)]
        z = sum(w)
        return [wi / z for wi in w]

    def mean_gap(lam):
        return _weighted_mean(x, posterior(lam)) - target

    # Bracket lambda: mean is increasing in lambda, so expand until it straddles.
    lo, hi = -1.0, 1.0
    while mean_gap(lo) > 0.0:
        lo *= 2.0
        if lo < -1e6:
            break
    while mean_gap(hi) < 0.0:
        hi *= 2.0
        if hi > 1e6:
            break

    lam = 0.0
    for _ in range(max_iter):
        lam = 0.5 * (lo + hi)
        g = mean_gap(lam)
        if abs(g) < tol:
            break
        if g > 0.0:
            hi = lam
        else:
            lo = lam
    return posterior(lam)


def relative_entropy(p, q):
    """Kullback-Leibler divergence ``sum_i p_i log(p_i / q_i)`` (nats).

    Non-negative, and zero exactly when ``p == q``. Terms with ``p_i == 0`` are
    dropped (limit ``0 log 0 = 0``).
    """
    if len(p) != len(q):
        raise ValueError("p and q must have the same length")
    total = 0.0
    for i in range(len(p)):
        if p[i] > 0.0:
            if q[i] <= 0.0:
                raise ValueError("q must be positive wherever p is positive")
            total += p[i] * math.log(p[i] / q[i])
    return total
