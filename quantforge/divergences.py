"""Divergences and distances between discrete probability distributions.

Ways to quantify how far apart two probability vectors ``p`` and ``q`` are:

  * ``kl_divergence`` -- ``sum p log(p/q)``: the information lost using ``q`` for ``p``.
    Asymmetric, unbounded, requires ``q_i > 0`` wherever ``p_i > 0``.
  * ``jensen_shannon_divergence`` -- symmetric, always finite KL to the mixture; its
    square root is a true metric, bounded by ``log 2``.
  * ``hellinger_distance`` -- ``(1/sqrt2) ||sqrt(p) - sqrt(q)||``: a bounded ``[0, 1]``
    metric.
  * ``total_variation_distance`` -- ``(1/2) sum |p - q|``: the largest probability gap
    over any event.
  * ``bhattacharyya_distance`` -- ``-log sum sqrt(p q)``: overlap-based, related to
    Hellinger.

Inputs are non-negative weight vectors; each is normalized to sum to one. Pure standard
library.
"""

import math


def _normalize(v):
    if any(x < 0 for x in v):
        raise ValueError("probabilities must be non-negative")
    s = sum(v)
    if s <= 0:
        raise ValueError("distribution must have positive total mass")
    return [x / s for x in v]


def _check(p, q):
    if len(p) != len(q):
        raise ValueError("p and q must have equal length")
    if len(p) == 0:
        raise ValueError("need at least one bin")
    return _normalize(p), _normalize(q)


def kl_divergence(p, q):
    """Kullback-Leibler divergence ``sum p_i log(p_i / q_i)`` in nats.

    Asymmetric and non-negative (zero iff ``p == q``). Raises if some ``q_i == 0`` where
    ``p_i > 0`` (the divergence is infinite there).
    """
    p, q = _check(p, q)
    total = 0.0
    for i in range(len(p)):
        if p[i] > 0.0:
            if q[i] <= 0.0:
                raise ValueError("q has zero mass where p is positive (KL is infinite)")
            total += p[i] * math.log(p[i] / q[i])
    return total


def jensen_shannon_divergence(p, q):
    """Jensen-Shannon divergence (symmetric, bounded by ``log 2``) in nats.

    ``0.5 KL(p || m) + 0.5 KL(q || m)`` with ``m = (p + q)/2``. Always finite; its
    square root is a metric.
    """
    p, q = _check(p, q)
    m = [(p[i] + q[i]) / 2.0 for i in range(len(p))]

    def _kl(a):
        t = 0.0
        for i in range(len(a)):
            if a[i] > 0.0:
                t += a[i] * math.log(a[i] / m[i])
        return t

    return 0.5 * _kl(p) + 0.5 * _kl(q)


def hellinger_distance(p, q):
    """Hellinger distance ``(1/sqrt2) sqrt(sum (sqrt(p_i) - sqrt(q_i))^2)`` in ``[0, 1]``."""
    p, q = _check(p, q)
    s = sum((math.sqrt(p[i]) - math.sqrt(q[i])) ** 2 for i in range(len(p)))
    return math.sqrt(s / 2.0)


def total_variation_distance(p, q):
    """Total-variation distance ``(1/2) sum |p_i - q_i|`` in ``[0, 1]``."""
    p, q = _check(p, q)
    return 0.5 * sum(abs(p[i] - q[i]) for i in range(len(p)))


def bhattacharyya_distance(p, q):
    """Bhattacharyya distance ``-log(sum sqrt(p_i q_i))`` (0 for identical, grows apart)."""
    p, q = _check(p, q)
    bc = sum(math.sqrt(p[i] * q[i]) for i in range(len(p)))
    if bc <= 0.0:
        return float("inf")
    return -math.log(bc)
