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
