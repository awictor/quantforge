"""Weighted descriptive statistics.

Summaries that give each observation a weight -- reliability weights, sampling weights,
or a probability distribution over scenarios. Provides the weighted mean, variance and
standard deviation (with the reliability-weight bias correction), and weighted quantiles
and median. Weights need not sum to 1; they are normalized internally. Pure standard
library.
"""


def weighted_mean(values, weights):
    """Weighted arithmetic mean ``sum(w x) / sum(w)``.

    ``weights`` must be non-negative and not all zero; they need not sum to 1.
    """
    if len(values) != len(weights):
        raise ValueError("values and weights must have equal length")
    if not values:
        raise ValueError("need at least one value")
    sw = 0.0
    swx = 0.0
    for x, w in zip(values, weights):
        if w < 0:
            raise ValueError("weights must be non-negative")
        sw += w
        swx += w * x
    if sw == 0.0:
        raise ValueError("weights sum to zero")
    return swx / sw


def weighted_variance(values, weights, unbiased=True):
    """Weighted variance about the weighted mean.

    With ``unbiased=True`` applies the reliability-weight correction
    ``V1 / (V1^2 - V2)`` where ``V1 = sum(w)`` and ``V2 = sum(w^2)`` (reduces to the
    ``1/(n-1)`` factor for equal weights); with ``unbiased=False`` divides by ``V1``
    (the population form). ``weights`` non-negative, not all zero.
    """
    if len(values) != len(weights):
        raise ValueError("values and weights must have equal length")
    if len(values) < 2:
        raise ValueError("need at least two values")
    mu = weighted_mean(values, weights)
    v1 = sum(weights)
    v2 = sum(w * w for w in weights)
    ss = sum(w * (x - mu) ** 2 for x, w in zip(values, weights))
    if unbiased:
        denom = v1 - v2 / v1
        if denom <= 0.0:
            raise ValueError("degenerate weights for unbiased variance")
        return ss / denom
    return ss / v1


def weighted_std(values, weights, unbiased=True):
    """Weighted standard deviation: square root of :func:`weighted_variance`."""
    return weighted_variance(values, weights, unbiased) ** 0.5


def weighted_quantile(values, weights, q):
    """Weighted ``q``-quantile (``q`` in ``[0, 1]``) by the cumulative-weight method.

    Sorts by value, forms the normalized cumulative weight at each point (midpoint
    convention), and linearly interpolates to the target ``q``. ``q=0.5`` is the weighted
    median. Reduces to the ordinary quantile when weights are equal.
    """
    if len(values) != len(weights):
        raise ValueError("values and weights must have equal length")
    if not values:
        raise ValueError("need at least one value")
    if not (0.0 <= q <= 1.0):
        raise ValueError("q must be in [0, 1]")
    pairs = sorted(zip(values, weights))
    total = sum(w for _, w in pairs)
    if total <= 0.0:
        raise ValueError("weights sum to zero")
    # Cumulative "plotting position" at each sorted point (midpoint of its weight band).
    cum = 0.0
    positions = []
    for x, w in pairs:
        if w < 0:
            raise ValueError("weights must be non-negative")
        positions.append((cum + w / 2.0) / total)
        cum += w
    xs = [x for x, _ in pairs]
    if q <= positions[0]:
        return xs[0]
    if q >= positions[-1]:
        return xs[-1]
    for i in range(1, len(positions)):
        if q <= positions[i]:
            p0, p1 = positions[i - 1], positions[i]
            frac = (q - p0) / (p1 - p0)
            return xs[i - 1] + frac * (xs[i] - xs[i - 1])
    return xs[-1]


def weighted_median(values, weights):
    """Weighted median: the weighted 0.5-quantile."""
    return weighted_quantile(values, weights, 0.5)
