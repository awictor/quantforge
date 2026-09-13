"""Backtests for Value-at-Risk and Expected Shortfall models.

Given a series of realized losses and the model's VaR (and ES) forecasts, these
tests check whether the risk model is well calibrated:

- ``kupiec_pof`` -- Kupiec's proportion-of-failures test: the exception rate should
  equal the VaR tail probability (unconditional coverage), a chi-square(1) LR test,
- ``christoffersen_independence`` -- exceptions should not cluster (independence), a
  chi-square(1) LR test on the exception indicator's Markov transitions,
- ``christoffersen_cc`` -- the joint conditional-coverage test (sum of the two, a
  chi-square(2)),
- ``acerbi_szekely_es`` -- the Acerbi-Szekely test statistic for Expected Shortfall
  calibration (near zero when ES is right, positive when ES understates the tail).

An exception is a period whose loss exceeds the VaR. Pure standard library on top
of the chi-square distribution.
"""

import math

from .distributions import chi2_sf


def _exceptions(losses, var_forecasts):
    n = len(losses)
    if n == 0 or len(var_forecasts) != n:
        raise ValueError("losses and var_forecasts must be equal-length, non-empty")
    return [1 if losses[i] > var_forecasts[i] else 0 for i in range(n)]


def kupiec_pof(losses, var_forecasts, alpha=0.01):
    """Kupiec proportion-of-failures (unconditional-coverage) test.

    ``alpha`` is the VaR tail probability (e.g. 0.01 for 99% VaR), so the expected
    exception rate is ``alpha``. Returns ``(LR, p_value)`` with ``LR`` a
    chi-square(1) likelihood ratio; a small p-value rejects the model's coverage.
    """
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be in (0, 1)")
    hits = _exceptions(losses, var_forecasts)
    n = len(hits)
    x = sum(hits)
    pi = x / n
    if x == 0 or x == n:
        # Degenerate MLE; fall back to the null term only.
        ll_null = x * math.log(alpha) + (n - x) * math.log(1.0 - alpha)
        lr = -2.0 * ll_null
        return lr, chi2_sf(lr, 1)
    ll_null = x * math.log(alpha) + (n - x) * math.log(1.0 - alpha)
    ll_alt = x * math.log(pi) + (n - x) * math.log(1.0 - pi)
    lr = -2.0 * (ll_null - ll_alt)
    return lr, chi2_sf(lr, 1)


def christoffersen_independence(losses, var_forecasts):
    """Christoffersen independence test: exceptions should not cluster.

    Fits a first-order Markov chain to the exception indicator and tests that the
    probability of an exception does not depend on whether the previous period was an
    exception. Returns ``(LR, p_value)``, a chi-square(1) LR test.
    """
    hits = _exceptions(losses, var_forecasts)
    n = len(hits)
    if n < 2:
        raise ValueError("need at least two observations")
    n00 = n01 = n10 = n11 = 0
    for t in range(1, n):
        prev, cur = hits[t - 1], hits[t]
        if prev == 0 and cur == 0:
            n00 += 1
        elif prev == 0 and cur == 1:
            n01 += 1
        elif prev == 1 and cur == 0:
            n10 += 1
        else:
            n11 += 1
    pi01 = n01 / (n00 + n01) if (n00 + n01) > 0 else 0.0
    pi11 = n11 / (n10 + n11) if (n10 + n11) > 0 else 0.0
    pi = (n01 + n11) / (n00 + n01 + n10 + n11)
    if pi <= 0.0 or pi >= 1.0:
        return 0.0, 1.0     # no exceptions or all exceptions: nothing to test

    def _ll(p):
        return 0.0 if p <= 0.0 or p >= 1.0 else None

    ll_null = ((n00 + n10) * math.log(1.0 - pi) + (n01 + n11) * math.log(pi))
    ll_alt = 0.0
    ll_alt += n00 * math.log(1.0 - pi01) if pi01 < 1.0 else 0.0
    ll_alt += n01 * math.log(pi01) if pi01 > 0.0 else 0.0
    ll_alt += n10 * math.log(1.0 - pi11) if pi11 < 1.0 else 0.0
    ll_alt += n11 * math.log(pi11) if pi11 > 0.0 else 0.0
    lr = -2.0 * (ll_null - ll_alt)
    lr = max(0.0, lr)
    return lr, chi2_sf(lr, 1)


def christoffersen_cc(losses, var_forecasts, alpha=0.01):
    """Christoffersen conditional-coverage test: coverage and independence jointly.

    The sum of :func:`kupiec_pof` and :func:`christoffersen_independence` LR
    statistics, tested as a chi-square(2). Returns ``(LR, p_value)``.
    """
    lr_uc, _ = kupiec_pof(losses, var_forecasts, alpha)
    lr_ind, _ = christoffersen_independence(losses, var_forecasts)
    lr = lr_uc + lr_ind
    return lr, chi2_sf(lr, 2)


def acerbi_szekely_es(losses, var_forecasts, es_forecasts, alpha=0.01):
    """Acerbi-Szekely (2014) Expected-Shortfall test statistic (their Z2).

    ``Z = (1 / (n alpha)) sum_t hit_t * loss_t / ES_t - 1``, where ``hit_t`` marks a
    VaR exception. Near zero when ES is well calibrated; positive when realized tail
    losses exceed the forecast ES (the model understates risk), negative when it
    overstates. Returns the scalar statistic.
    """
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be in (0, 1)")
    n = len(losses)
    if n == 0 or len(var_forecasts) != n or len(es_forecasts) != n:
        raise ValueError("losses, var_forecasts, es_forecasts must align, non-empty")
    if any(e <= 0.0 for e in es_forecasts):
        raise ValueError("ES forecasts must be positive")
    hits = _exceptions(losses, var_forecasts)
    s = sum(hits[t] * losses[t] / es_forecasts[t] for t in range(n))
    return s / (n * alpha) - 1.0
