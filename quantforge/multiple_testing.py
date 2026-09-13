"""Multiple-hypothesis-testing corrections.

Running many tests inflates the chance of a false positive. These routines adjust a
set of raw p-values so a fixed threshold controls a family-wise error rate or a
false-discovery rate:

- ``bonferroni`` -- multiply by ``m`` (controls the family-wise error rate; simplest
  and most conservative),
- ``holm`` -- step-down family-wise control, uniformly more powerful than
  Bonferroni,
- ``benjamini_hochberg`` -- step-up control of the false-discovery rate under
  independence or positive dependence,
- ``benjamini_yekutieli`` -- FDR control under arbitrary dependence (adds the
  harmonic-number penalty).

Each returns a list of adjusted p-values aligned with the input, clamped to
``[0, 1]`` and enforced monotone so a fixed cutoff gives coherent rejections. Pure
standard library.
"""


def _check(pvals):
    if not pvals:
        raise ValueError("need at least one p-value")
    if any(not (0.0 <= p <= 1.0) for p in pvals):
        raise ValueError("p-values must be in [0, 1]")


def bonferroni(pvals):
    """Bonferroni-adjusted p-values: ``min(1, m * p)`` for ``m`` tests."""
    _check(pvals)
    m = len(pvals)
    return [min(1.0, m * p) for p in pvals]


def holm(pvals):
    """Holm step-down family-wise adjusted p-values.

    Sorts ascending, scales the ``k``-th smallest (0-based) by ``m - k``, then takes
    a running maximum so the sequence is monotone, and unshuffles to the input
    order. Controls the family-wise error rate and dominates :func:`bonferroni`.
    """
    _check(pvals)
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running = 0.0
    for rank, idx in enumerate(order):
        val = min(1.0, (m - rank) * pvals[idx])
        running = max(running, val)
        adj[idx] = running
    return adj


def benjamini_hochberg(pvals):
    """Benjamini-Hochberg FDR-adjusted p-values (step-up).

    Sorts ascending, scales the ``k``-th smallest (1-based) by ``m / k``, takes a
    running minimum from the largest down so the sequence is monotone, and unshuffles
    to the input order. Controls the false-discovery rate under independence or
    positive dependence; less conservative than family-wise methods.
    """
    _check(pvals)
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running = 1.0
    for rank in range(m - 1, -1, -1):
        idx = order[rank]
        val = min(1.0, pvals[idx] * m / (rank + 1))
        running = min(running, val)
        adj[idx] = running
    return adj


def benjamini_yekutieli(pvals):
    """Benjamini-Yekutieli FDR-adjusted p-values (arbitrary dependence).

    As Benjamini-Hochberg but with the extra factor ``c(m) = sum_{i=1}^m 1/i``, which
    makes the procedure valid under any dependence structure at the cost of power.
    """
    _check(pvals)
    m = len(pvals)
    c = sum(1.0 / i for i in range(1, m + 1))
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running = 1.0
    for rank in range(m - 1, -1, -1):
        idx = order[rank]
        val = min(1.0, pvals[idx] * m * c / (rank + 1))
        running = min(running, val)
        adj[idx] = running
    return adj
