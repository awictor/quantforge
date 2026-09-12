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


def panjer_negative_binomial(size, prob, severity_pmf, max_k=None):
    """Aggregate-loss distribution for compound negative-binomial claim counts.

    The claim count ``N ~ NegBinom(size=r, prob=p)`` (``P(N=n) = C(n+r-1, n)
    p^r (1-p)^n``, mean ``r(1-p)/p``) is over-dispersed relative to Poisson
    (variance > mean), capturing claim contagion. It is a Panjer ``(a, b)`` class
    with ``a = 1 - p`` and ``b = (r - 1)(1 - p)``:

        g_0 = p^r  (if severity has no mass at 0),
        g_k = 1/(1 - a f_0) * sum_{j=1}^{k} (a + b j / k) f_j g_{k-j}.

    Parameters
    ----------
    size : float
        The NB ``r`` (number of failures); ``r > 0``.
    prob : float
        The NB success probability ``p`` in ``(0, 1]``.
    severity_pmf : sequence of float
        Severity probabilities on an integer grid.
    max_k : int, optional
        Aggregate grid cutoff.

    Returns
    -------
    list[float]
        ``g[k] = P(S = k)``; sums to ~1.
    """
    import math
    if size <= 0:
        raise ValueError("size must be positive")
    if not (0.0 < prob <= 1.0):
        raise ValueError("prob must be in (0, 1]")
    n = len(severity_pmf)
    if n == 0:
        raise ValueError("severity_pmf must be non-empty")
    s = sum(severity_pmf)
    if s <= 0:
        raise ValueError("severity_pmf must sum to a positive value")
    f = [p / s for p in severity_pmf]
    a = 1.0 - prob
    b = (size - 1.0) * (1.0 - prob)
    mean_n = size * (1.0 - prob) / prob
    if max_k is None:
        max_k = int(mean_n * n) * 4 + 20

    g = [0.0] * (max_k + 1)
    # P(S=0) = P(N=0) when severity has no atom at 0; general PGF at f_0.
    g[0] = prob ** size / (1.0 - a * f[0]) ** size
    denom = 1.0 - a * f[0]
    for k in range(1, max_k + 1):
        acc = 0.0
        for j in range(1, min(k, n - 1) + 1):
            acc += (a + b * j / k) * f[j] * g[k - j]
        g[k] = acc / denom
    return g
