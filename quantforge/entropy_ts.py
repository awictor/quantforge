"""Time-series regularity measures: approximate, sample and permutation entropy.

These quantify how unpredictable a series is:

- ``approximate_entropy`` (Pincus) -- the log-likelihood that runs close for ``m``
  points stay close for ``m + 1``; low for regular signals, high for irregular ones,
  but biased by self-matches on short series,
- ``sample_entropy`` (Richman-Moorman) -- the same idea without self-matches, so it
  is largely length-independent and the usual choice,
- ``permutation_entropy`` (Bandt-Pompe) -- the Shannon entropy of the ordinal
  patterns of length ``m``, normalized to ``[0, 1]``; robust to monotone transforms
  and to noise.

Pure standard library.
"""

import math


def _phi(series, m, r):
    """Fraction-of-close-matches term used by approximate entropy."""
    n = len(series)
    patterns = [series[i:i + m] for i in range(n - m + 1)]
    total = 0.0
    for i in range(len(patterns)):
        count = 0
        for j in range(len(patterns)):
            if max(abs(patterns[i][k] - patterns[j][k]) for k in range(m)) <= r:
                count += 1
        total += math.log(count / (n - m + 1))
    return total / (n - m + 1)


def approximate_entropy(series, m=2, r=None):
    """Approximate entropy ``ApEn(m, r)`` (Pincus).

    ``r`` is the tolerance for a match (defaults to ``0.2 * std(series)``). Returns
    ``phi_m - phi_{m+1}``; near zero for a perfectly regular series and larger for an
    irregular one. Includes self-matches, so it is biased low on short series -- use
    :func:`sample_entropy` when that matters.
    """
    n = len(series)
    if n < m + 2:
        raise ValueError("series too short for the chosen m")
    if r is None:
        mean = sum(series) / n
        sd = math.sqrt(sum((x - mean) ** 2 for x in series) / n)
        r = 0.2 * sd
    if r <= 0.0:
        raise ValueError("tolerance r must be positive (constant series?)")
    return _phi(series, m, r) - _phi(series, m + 1, r)


def sample_entropy(series, m=2, r=None):
    """Sample entropy ``SampEn(m, r)`` (Richman-Moorman).

    ``-log(A / B)`` where ``B`` counts template matches of length ``m`` and ``A`` of
    length ``m + 1``, both excluding self-matches. Larger means less regular / more
    complex. ``r`` defaults to ``0.2 * std(series)``. Raises if no longer-template
    matches occur (entropy would be infinite).
    """
    n = len(series)
    if n < m + 2:
        raise ValueError("series too short for the chosen m")
    if r is None:
        mean = sum(series) / n
        sd = math.sqrt(sum((x - mean) ** 2 for x in series) / n)
        r = 0.2 * sd
    if r <= 0.0:
        raise ValueError("tolerance r must be positive (constant series?)")

    def count(mm):
        pats = [series[i:i + mm] for i in range(n - mm + 1)]
        c = 0
        for i in range(len(pats)):
            for j in range(i + 1, len(pats)):
                if max(abs(pats[i][k] - pats[j][k]) for k in range(mm)) <= r:
                    c += 1
        return c

    b = count(m)
    a = count(m + 1)
    if b == 0 or a == 0:
        raise ValueError("no template matches; increase r or the series length")
    return -math.log(a / b)


def permutation_entropy(series, m=3, normalize=True):
    """Permutation entropy (Bandt-Pompe) from the ordinal patterns of length ``m``.

    Slides a window of ``m`` points, records which permutation sorts each window, and
    returns the Shannon entropy of the pattern distribution. With ``normalize=True``
    it is divided by ``log(m!)`` to land in ``[0, 1]`` -- 0 for a monotone series, 1
    for one whose orderings are uniformly random. Invariant to any monotone transform
    of the series.
    """
    n = len(series)
    if n < m + 1:
        raise ValueError("series too short for the chosen m")
    if m < 2:
        raise ValueError("m must be at least 2")
    counts = {}
    for i in range(n - m + 1):
        window = series[i:i + m]
        # Ordinal pattern: the permutation of indices that sorts the window.
        pattern = tuple(sorted(range(m), key=lambda k: window[k]))
        counts[pattern] = counts.get(pattern, 0) + 1
    total = sum(counts.values())
    h = 0.0
    for c in counts.values():
        p = c / total
        h -= p * math.log(p)
    if normalize:
        h /= math.log(math.factorial(m))
    return h
