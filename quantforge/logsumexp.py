"""Numerically stable log-sum-exp, softmax, and log-softmax.

``log(sum(exp(x)))`` overflows the moment any ``x`` is large and underflows to
``-inf`` when all are very negative. The standard fix subtracts the maximum first,

    logsumexp(x) = m + log(sum(exp(x - m))),   m = max(x),

which is exact and never overflows. ``softmax`` and ``log_softmax`` build on it -- the
normalized-probability and log-probability transforms at the heart of classification,
mixture weights, and Boltzmann-style models. Optional ``weights`` give a
log-weighted-sum-exp. Pure standard library.
"""

import math


def logsumexp(x, weights=None):
    """Stable ``log(sum_i w_i exp(x_i))`` (weights default to 1).

    Shifts by the maximum so it never overflows; returns ``-inf`` if every weighted
    term is zero. With ``weights`` it is the log of a weighted sum of exponentials
    (weights must be non-negative).
    """
    if not x:
        raise ValueError("need at least one value")
    if weights is not None:
        if len(weights) != len(x):
            raise ValueError("weights must match x in length")
        if any(w < 0 for w in weights):
            raise ValueError("weights must be non-negative")
    m = max(x)
    if m == float("-inf"):
        return float("-inf")
    if weights is None:
        s = sum(math.exp(xi - m) for xi in x)
    else:
        s = sum(w * math.exp(xi - m) for xi, w in zip(x, weights))
    if s <= 0.0:
        return float("-inf")
    return m + math.log(s)


def softmax(x):
    """Stable softmax: ``exp(x_i) / sum_j exp(x_j)``, a probability vector summing to 1."""
    if not x:
        raise ValueError("need at least one value")
    m = max(x)
    exps = [math.exp(xi - m) for xi in x]
    total = sum(exps)
    return [e / total for e in exps]


def log_softmax(x):
    """Stable log-softmax: ``x_i - logsumexp(x)`` (avoids the overflow of ``log(softmax)``)."""
    if not x:
        raise ValueError("need at least one value")
    lse = logsumexp(x)
    return [xi - lse for xi in x]
