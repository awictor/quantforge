"""Panjer recursion for compound-Poisson aggregate loss distributions.

The aggregate loss ``S = X_1 + ... + X_N`` of ``N ~ Poisson(lam)`` claims, each an
i.i.d. severity on a discrete grid, has a distribution computed exactly by Panjer's
recursion (no simulation). For Poisson frequency (``a = 0``, ``b = lam``):

    g_0 = exp(-lam (1 - f_0)),
    g_k = (lam / k) * sum_{j=1}^{k} j f_j g_{k-j},

where ``f_j`` is the severity probability at grid point ``j`` and ``g_k`` the
aggregate probability. From the aggregate distribution the mean, stop-loss
premiums, and layer expected losses follow directly. Pure standard library.
"""


def panjer_poisson(lam, severity_pmf, max_k=None):
    """Aggregate-loss distribution for compound Poisson via Panjer recursion.

    Parameters
    ----------
    lam : float
        Poisson claim frequency (mean number of claims).
    severity_pmf : sequence of float
        Severity probabilities on an integer grid ``0, 1, 2, ...`` (index =
        severity in grid units); should sum to 1.
    max_k : int, optional
        Highest aggregate grid point to compute. Defaults to a cutoff capturing
        essentially all mass (``ceil(lam * n) * 4 + 20``).

    Returns
    -------
    list[float]
        ``g[k] = P(S = k)`` on the aggregate grid; sums to ~1.
    """
    if lam < 0:
        raise ValueError("lam must be non-negative")
    n = len(severity_pmf)
    if n == 0:
        raise ValueError("severity_pmf must be non-empty")
    s = sum(severity_pmf)
    if s <= 0:
        raise ValueError("severity_pmf must sum to a positive value")
    f = [p / s for p in severity_pmf]
    if max_k is None:
        max_k = int(lam * n) * 4 + 20

    import math
    g = [0.0] * (max_k + 1)
    g[0] = math.exp(-lam * (1.0 - f[0]))
    for k in range(1, max_k + 1):
        acc = 0.0
        for j in range(1, min(k, n - 1) + 1):
            acc += j * f[j] * g[k - j]
        g[k] = (lam / k) * acc
    return g


def aggregate_mean(g):
    """Mean of an aggregate distribution ``g`` (grid units)."""
    return sum(k * g[k] for k in range(len(g)))


def stop_loss_premium(g, retention):
    """Stop-loss premium ``E[max(S - retention, 0)]`` from the aggregate grid."""
    return sum((k - retention) * g[k] for k in range(len(g)) if k > retention)


def layer_expected_loss(g, attachment, limit):
    """Expected loss to a reinsurance layer ``[attachment, attachment + limit]``.

    ``E[min(max(S - attachment, 0), limit)]`` -- the excess-of-loss layer cost.
    """
    top = attachment + limit
    total = 0.0
    for k in range(len(g)):
        loss = min(max(k - attachment, 0.0), limit)
        total += loss * g[k]
    return total
