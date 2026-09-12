"""Buhlmann and Buhlmann-Straub credibility (experience rating).

Credibility blends a risk's own experience with the collective mean:

    premium = Z * own_mean + (1 - Z) * collective_mean,

where the credibility factor ``Z = n / (n + k)`` weights the individual data by
its volume ``n`` against a stiffness ``k = EPV / VHM`` -- the ratio of the
expected process variance (within-risk noise) to the variance of hypothetical
means (between-risk spread). More data, or a larger between-risk spread relative to
noise, pulls ``Z`` toward 1 (trust the individual); more noise pulls it toward 0
(trust the collective). Pure standard library.
"""


def buhlmann_k(expected_process_variance, variance_of_hypothetical_means):
    """Buhlmann stiffness ``k = EPV / VHM``.

    ``EPV`` is the mean within-risk variance; ``VHM`` the between-risk variance of
    the true means. Smaller ``k`` -> more credibility to the individual.
    """
    if variance_of_hypothetical_means <= 0.0:
        raise ValueError("VHM must be positive")
    if expected_process_variance < 0.0:
        raise ValueError("EPV must be non-negative")
    return expected_process_variance / variance_of_hypothetical_means


def credibility_factor(n, k):
    """Buhlmann credibility ``Z = n / (n + k)`` for ``n`` observations."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if k < 0:
        raise ValueError("k must be non-negative")
    if n + k == 0:
        return 0.0
    return n / (n + k)


def buhlmann_premium(own_mean, collective_mean, n, epv, vhm):
    """Buhlmann credibility premium.

    Blends the risk's own mean with the collective mean using ``Z = n/(n+k)``,
    ``k = EPV/VHM``. Lies between the two means; approaches the own mean as ``n``
    grows or the between-risk spread dominates the noise.
    """
    k = buhlmann_k(epv, vhm)
    z = credibility_factor(n, k)
    return z * own_mean + (1.0 - z) * collective_mean


def buhlmann_straub_premium(claims, exposures, collective_mean, epv, vhm):
    """Buhlmann-Straub premium for a risk with per-period exposures.

    Generalizes Buhlmann to unequal exposures ``m_i``: the credibility uses total
    exposure ``m = sum m_i`` with ``Z = m / (m + k)`` and the own estimate is the
    exposure-weighted claim rate ``sum claims_i / m``.
    """
    if len(claims) != len(exposures):
        raise ValueError("claims and exposures must have the same length")
    m = sum(exposures)
    if m <= 0.0:
        raise ValueError("total exposure must be positive")
    own_rate = sum(claims) / m
    k = buhlmann_k(epv, vhm)
    z = m / (m + k)
    return z * own_rate + (1.0 - z) * collective_mean
