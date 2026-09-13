"""Wald-Wolfowitz runs test for randomness of a sequence.

A *run* is a maximal streak of like values. Too few runs means the sequence
clusters (trending / positive dependence); too many means it over-alternates
(mean-reverting). Under randomness the run count is approximately normal with a
known mean and variance, giving a two-sided z-test.

``runs_test`` dichotomizes a numeric series about its median (values equal to the
median are dropped); ``runs_test_binary`` takes a 0/1 (or two-symbol) sequence
directly. Both return ``(z, p_value)``. Pure standard library.
"""

import math

from .mathfns import norm_cdf


def _runs_z(n1, n2, runs):
    """Two-sided z-statistic and p-value for ``runs`` given group sizes n1, n2."""
    n = n1 + n2
    if n1 == 0 or n2 == 0:
        raise ValueError("need both symbols present")
    mean = 2.0 * n1 * n2 / n + 1.0
    var = (2.0 * n1 * n2 * (2.0 * n1 * n2 - n)) / (n * n * (n - 1.0))
    if var <= 0.0:
        return 0.0, 1.0
    z = (runs - mean) / math.sqrt(var)
    p = 2.0 * (1.0 - norm_cdf(abs(z)))
    return z, p


def _count_runs(symbols):
    runs = 1
    for i in range(1, len(symbols)):
        if symbols[i] != symbols[i - 1]:
            runs += 1
    return runs


def runs_test_binary(sequence):
    """Wald-Wolfowitz runs test on a two-symbol sequence.

    ``sequence`` is any list of two distinct values (e.g. 0/1, +/-). Returns
    ``(z, p_value)``; a negative ``z`` (few runs) signals clustering, a positive
    ``z`` (many runs) over-alternation. Requires at least one of each symbol.
    """
    if len(sequence) < 2:
        raise ValueError("need at least two elements")
    distinct = list(dict.fromkeys(sequence))
    if len(distinct) != 2:
        raise ValueError("sequence must contain exactly two distinct symbols")
    n1 = sum(1 for s in sequence if s == distinct[0])
    n2 = len(sequence) - n1
    return _runs_z(n1, n2, _count_runs(sequence))


def runs_test(values):
    """Runs test on a numeric series, dichotomized about its median.

    Values above the median are one symbol, below the other; values exactly equal to
    the median are dropped. Returns ``(z, p_value)`` for the null that the sequence
    of above/below signs is random. A small p-value rejects randomness: ``z < 0`` for
    trending/clustered data, ``z > 0`` for over-alternating (mean-reverting) data.
    """
    n = len(values)
    if n < 2:
        raise ValueError("need at least two values")
    s = sorted(values)
    m = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
    signs = [1 if v > m else 0 for v in values if v != m]
    if len(signs) < 2:
        raise ValueError("too many values equal the median")
    n1 = sum(signs)
    n2 = len(signs) - n1
    return _runs_z(n1, n2, _count_runs(signs))
