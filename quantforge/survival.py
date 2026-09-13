"""Nonparametric survival estimators for right-censored data.

Given event times, some of which are *censored* (the subject left the study or a
loan matured before default), these estimate the survival and cumulative-hazard
curves without assuming a parametric form:

- ``kaplan_meier`` -- the product-limit estimator of the survival function
  ``S(t) = prod_{t_i <= t} (1 - d_i / n_i)``, where at each event time ``d_i`` events
  occur among ``n_i`` still at risk,
- ``nelson_aalen`` -- the cumulative hazard ``H(t) = sum_{t_i <= t} d_i / n_i``.

The two are linked by ``S(t) ~ exp(-H(t))``. Pure standard library.
"""

import math


def _risk_table(times, events):
    """Distinct event times with (d_i events, n_i at risk), ascending.

    ``events[i] = 1`` marks an observed event, ``0`` a right-censored observation.
    Censored subjects leave the risk set but contribute no event.
    """
    n = len(times)
    if n == 0 or len(events) != n:
        raise ValueError("times and events must be equal-length, non-empty")
    if any(t < 0 for t in times):
        raise ValueError("times must be non-negative")
    order = sorted(range(n), key=lambda i: times[i])
    rows = []
    i = 0
    at_risk = n
    while i < n:
        t = times[order[i]]
        d = 0
        censored = 0
        j = i
        while j < n and times[order[j]] == t:
            if events[order[j]]:
                d += 1
            else:
                censored += 1
            j += 1
        rows.append((t, d, censored, at_risk))
        at_risk -= (d + censored)
        i = j
    return rows


def kaplan_meier(times, events):
    """Kaplan-Meier product-limit survival estimate.

    Returns ``(event_times, survival)`` giving the step-function survival ``S(t)``
    at each distinct time where an event occurs; censoring times only shrink the risk
    set. ``S`` starts at 1, is non-increasing, and drops by the factor
    ``1 - d_i / n_i`` at each event time. With no censoring it equals
    ``1 - ECDF(t)``.
    """
    rows = _risk_table(times, events)
    out_t = []
    out_s = []
    s = 1.0
    for t, d, _c, n_i in rows:
        if d == 0:
            continue                     # pure censoring time: no jump
        s *= (1.0 - d / n_i)
        out_t.append(t)
        out_s.append(s)
    return out_t, out_s


def nelson_aalen(times, events):
    """Nelson-Aalen cumulative-hazard estimate.

    Returns ``(event_times, cumulative_hazard)`` with ``H(t) = sum d_i / n_i`` over
    event times up to ``t``. Non-decreasing from 0; ``exp(-H(t))`` approximates the
    Kaplan-Meier survival (they agree closely when the per-step hazard is small).
    """
    rows = _risk_table(times, events)
    out_t = []
    out_h = []
    h = 0.0
    for t, d, _c, n_i in rows:
        if d == 0:
            continue
        h += d / n_i
        out_t.append(t)
        out_h.append(h)
    return out_t, out_h


def survival_at(times, events, query, estimator="km"):
    """Evaluate the survival function at ``query`` from a fitted step curve.

    ``estimator`` is ``"km"`` (Kaplan-Meier) or ``"na"`` (``exp(-Nelson-Aalen)``).
    Returns the survival at the largest event time ``<= query`` (1 before the first
    event).
    """
    if estimator == "km":
        ts, ss = kaplan_meier(times, events)
    elif estimator == "na":
        ts, hs = nelson_aalen(times, events)
        ss = [math.exp(-h) for h in hs]
    else:
        raise ValueError("estimator must be 'km' or 'na'")
    s = 1.0
    for i in range(len(ts)):
        if ts[i] <= query:
            s = ss[i]
        else:
            break
    return s


def log_rank_test(times1, events1, times2, events2):
    """Log-rank (Mantel-Cox) test comparing two survival curves.

    At each distinct event time across the pooled sample, compares the observed
    events in group 1 with the number expected under the null of equal hazards
    (proportional to each group's share of the risk set), accumulating the
    observed-minus-expected and its hypergeometric variance. The statistic

        chi2 = (sum (O1 - E1))^2 / sum V1

    is asymptotically chi-square(1). Returns ``(chi2, p_value)``; a small p-value
    rejects equal survival between the groups. Uses the chi-square survival function.
    """
    from .distributions import chi2_sf

    rows1 = {t: (d, c, n) for t, d, c, n in _risk_table(times1, events1)}
    rows2 = {t: (d, c, n) for t, d, c, n in _risk_table(times2, events2)}
    all_times = sorted(set(rows1) | set(rows2))

    n1 = len(times1)
    n2 = len(times2)
    if n1 == 0 or n2 == 0:
        raise ValueError("both groups must be non-empty")

    # Walk the pooled timeline, tracking each group's risk set.
    risk1, risk2 = n1, n2
    # Map time -> (events, censored) for each group.
    ec1 = {t: (rows1[t][0], rows1[t][1]) for t in rows1}
    ec2 = {t: (rows2[t][0], rows2[t][1]) for t in rows2}

    o_minus_e = 0.0
    var = 0.0
    for t in all_times:
        d1, c1 = ec1.get(t, (0, 0))
        d2, c2 = ec2.get(t, (0, 0))
        d = d1 + d2
        n = risk1 + risk2
        if n > 1 and d > 0:
            e1 = d * risk1 / n
            o_minus_e += d1 - e1
            var += (d * (risk1 / n) * (1 - risk1 / n) * (n - d) / (n - 1))
        risk1 -= (d1 + c1)
        risk2 -= (d2 + c2)

    if var <= 0.0:
        return 0.0, 1.0
    chi2 = o_minus_e * o_minus_e / var
    return chi2, chi2_sf(chi2, 1)
