"""Benford's law first-digit analysis for anomaly / fraud detection.

Naturally-occurring numbers that span several orders of magnitude have leading
digits distributed as ``P(d) = log10(1 + 1/d)`` -- digit 1 leads ~30.1% of the time,
digit 9 only ~4.6%. Fabricated or constrained data often departs from this, so a
conformance test flags datasets worth a second look:

- ``benford_expected`` -- the reference first-digit probabilities,
- ``first_digit`` / ``first_digit_distribution`` -- extract and tabulate leading
  digits,
- ``benford_chi_square`` -- chi-square goodness-of-fit against Benford (8 df),
- ``benford_mad`` -- Nigrini's mean absolute deviation conformance statistic.

Pure standard library on top of the chi-square distribution.
"""

import math


def benford_expected():
    """First-digit probabilities under Benford's law, ``d = 1..9``.

    Returns a list of 9 probabilities ``log10(1 + 1/d)`` summing to one.
    """
    return [math.log10(1.0 + 1.0 / d) for d in range(1, 10)]


def first_digit(x):
    """Leading (most significant) decimal digit of ``x``, ignoring sign and zeros.

    Returns an integer ``1..9``; raises for zero (no leading digit).
    """
    x = abs(float(x))
    if x == 0.0:
        raise ValueError("zero has no leading digit")
    # Scale into [1, 10).
    while x >= 10.0:
        x /= 10.0
    while x < 1.0:
        x *= 10.0
    return int(x)


def first_digit_distribution(values):
    """Observed first-digit counts and proportions for a dataset.

    Returns ``(counts, proportions)``, each a list of length 9 for digits ``1..9``.
    Zero values are skipped. Requires at least one non-zero value.
    """
    counts = [0] * 9
    total = 0
    for v in values:
        if v == 0:
            continue
        counts[first_digit(v) - 1] += 1
        total += 1
    if total == 0:
        raise ValueError("no non-zero values")
    proportions = [c / total for c in counts]
    return counts, proportions


def benford_chi_square(values):
    """Chi-square goodness-of-fit of first digits against Benford's law.

    ``sum (O_d - E_d)^2 / E_d`` over digits ``1..9`` with ``E_d = n P_benford(d)``,
    referenced to a chi-square with 8 degrees of freedom. Returns
    ``(statistic, p_value)``; a small p-value rejects Benford conformance.
    """
    from .distributions import chi2_sf

    counts, _ = first_digit_distribution(values)
    n = sum(counts)
    exp = benford_expected()
    stat = sum((counts[i] - n * exp[i]) ** 2 / (n * exp[i]) for i in range(9))
    return stat, chi2_sf(stat, 8)


def benford_mad(values):
    """Nigrini's mean absolute deviation from Benford's first-digit law.

    ``MAD = (1/9) sum_d |observed_prop(d) - benford_prop(d)|``. Nigrini's rule of
    thumb: below ~0.006 is close conformance, above ~0.015 is nonconformity.
    Independent of sample size, unlike the chi-square statistic.
    """
    _, prop = first_digit_distribution(values)
    exp = benford_expected()
    return sum(abs(prop[i] - exp[i]) for i in range(9)) / 9.0
