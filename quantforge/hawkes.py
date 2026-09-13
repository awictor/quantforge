"""Self-exciting (Hawkes) point process with an exponential kernel.

A Hawkes process models event clustering -- trades, order arrivals, defaults -- where
each event raises the arrival intensity of the next. With an exponential kernel the
conditional intensity is

    lambda(t) = mu + sum_{t_i < t} alpha * exp(-beta (t - t_i)),

with baseline ``mu > 0``, jump ``alpha >= 0`` and decay ``beta > alpha``. The ratio
``n = alpha / beta`` is the branching ratio (the expected number of children per
event); ``n < 1`` is required for stationarity. This module gives the intensity, the
exact recursive log-likelihood, Ogata thinning simulation, and a maximum-likelihood
fit. Pure standard library.
"""

import math

from .optimize import nelder_mead


def intensity(t, history, mu, alpha, beta):
    """Conditional intensity ``lambda(t)`` given past event times ``history``.

    Sums the exponential kernel over events strictly before ``t``.
    """
    if mu <= 0.0 or alpha < 0.0 or beta <= 0.0:
        raise ValueError("require mu > 0, alpha >= 0, beta > 0")
    s = 0.0
    for ti in history:
        if ti < t:
            s += math.exp(-beta * (t - ti))
    return mu + alpha * s


def branching_ratio(alpha, beta):
    """Branching ratio ``alpha / beta`` -- expected offspring per event.

    Below 1 the process is stationary; at or above 1 it explodes.
    """
    if beta <= 0.0:
        raise ValueError("beta must be positive")
    return alpha / beta


def log_likelihood(events, mu, alpha, beta, t_end=None):
    """Exact log-likelihood of an exponential-kernel Hawkes process.

    Uses the standard recursion for the excitation term ``A_i = sum_{j<i}
    exp(-beta (t_i - t_j))`` (updated as ``A_i = exp(-beta dt)(1 + A_{i-1})``), so it
    is linear in the number of events. ``events`` is a sorted list of event times on
    ``[0, t_end]`` (``t_end`` defaults to the last event). The log-likelihood is

        sum_i log(mu + alpha A_i) - mu T - (alpha/beta) sum_i (1 - exp(-beta (T-t_i))).
    """
    if mu <= 0.0 or alpha < 0.0 or beta <= 0.0:
        raise ValueError("require mu > 0, alpha >= 0, beta > 0")
    n = len(events)
    if n == 0:
        raise ValueError("need at least one event")
    T = events[-1] if t_end is None else t_end
    ll = 0.0
    a = 0.0
    prev = events[0]
    for i in range(n):
        if i == 0:
            a = 0.0
        else:
            a = math.exp(-beta * (events[i] - prev)) * (1.0 + a)
            prev = events[i]
        ll += math.log(mu + alpha * a)
    # Compensator (integral of the intensity over [0, T]).
    comp = mu * T
    for ti in events:
        comp += (alpha / beta) * (1.0 - math.exp(-beta * (T - ti)))
    return ll - comp


def simulate(mu, alpha, beta, t_max, seed=1234567):
    """Simulate a Hawkes process on ``[0, t_max]`` by Ogata's thinning algorithm.

    Returns the sorted list of event times. Requires ``alpha < beta`` (stationarity)
    for a well-behaved simulation. Deterministic given ``seed`` (an LCG uniform
    stream).
    """
    if mu <= 0.0 or alpha < 0.0 or beta <= 0.0:
        raise ValueError("require mu > 0, alpha >= 0, beta > 0")
    state = seed & 0x7FFFFFFF

    def rand():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return (state + 0.5) / 0x80000000

    events = []
    t = 0.0
    while t < t_max:
        # Upper bound: intensity just after t, including any event exactly at t
        # (the kernel decays between events, so this bounds the next interval).
        lam_bar = mu + alpha * sum(math.exp(-beta * (t - ti))
                                   for ti in events if ti <= t)
        u = rand()
        w = -math.log(u) / lam_bar          # candidate inter-arrival at upper bound
        t += w
        if t >= t_max:
            break
        d = rand()
        if d * lam_bar <= intensity(t, events, mu, alpha, beta):
            events.append(t)
    return events


def fit(events, t_end=None, x0=None):
    """Maximum-likelihood fit of ``(mu, alpha, beta)`` to observed event times.

    Maximizes :func:`log_likelihood` with Nelder-Mead over a log-parameterization
    (keeping ``mu, alpha, beta`` positive). Returns a dict with ``mu``, ``alpha``,
    ``beta``, ``branching_ratio`` and the attained ``log_likelihood``. Provide
    ``x0 = (mu, alpha, beta)`` to seed the search.
    """
    n = len(events)
    if n < 2:
        raise ValueError("need at least two events")
    T = events[-1] if t_end is None else t_end
    if x0 is None:
        base = n / T
        x0 = (0.5 * base, 0.5 * base, base)

    def neg_ll(theta):
        mu, alpha, beta = (math.exp(theta[0]), math.exp(theta[1]), math.exp(theta[2]))
        if alpha >= beta:            # enforce stationarity
            return 1e12
        try:
            return -log_likelihood(events, mu, alpha, beta, T)
        except (ValueError, OverflowError):
            return 1e12

    start = [math.log(x0[0]), math.log(x0[1]), math.log(x0[2])]
    best_x, _ = nelder_mead(neg_ll, start, max_iter=4000)
    mu, alpha, beta = math.exp(best_x[0]), math.exp(best_x[1]), math.exp(best_x[2])
    return {
        "mu": mu,
        "alpha": alpha,
        "beta": beta,
        "branching_ratio": alpha / beta,
        "log_likelihood": log_likelihood(events, mu, alpha, beta, T),
    }


def residuals(events, mu, alpha, beta):
    """Time-rescaling residuals of a fitted Hawkes process.

    By the time-rescaling theorem the integrated intensity (compensator) between
    consecutive events, ``tau_i = integral_{t_{i-1}}^{t_i} lambda(u) du``, is a
    sequence of i.i.d. unit-rate exponentials when the model is correctly specified.
    Computes those ``tau_i`` for ``i >= 1`` using the exponential-kernel recursion,

        tau_i = mu (t_i - t_{i-1})
                + (alpha/beta) sum_{j<i} [exp(-beta (t_{i-1}-t_j)) - exp(-beta (t_i-t_j))],

    accumulated via a running kernel sum. Returns the list of ``n - 1`` residuals.
    """
    if mu <= 0.0 or alpha < 0.0 or beta <= 0.0:
        raise ValueError("require mu > 0, alpha >= 0, beta > 0")
    n = len(events)
    if n < 2:
        raise ValueError("need at least two events")
    taus = []
    # B_i = sum_{j : t_j <= t_{i-1}} exp(-beta (t_{i-1} - t_j)), the kernel sum at the
    # left endpoint of interval i (including the event at t_{i-1}). B_1 = 1 (only t_0),
    # and B_{i+1} = 1 + exp(-beta dt_i) B_i.
    b = 1.0
    for i in range(1, n):
        dt = events[i] - events[i - 1]
        tau = mu * dt + (alpha / beta) * (1.0 - math.exp(-beta * dt)) * b
        taus.append(tau)
        b = 1.0 + math.exp(-beta * dt) * b
    return taus


def gof_test(events, mu, alpha, beta):
    """Kolmogorov-Smirnov goodness-of-fit of a fitted Hawkes model.

    Applies the time-rescaling theorem: under a correct model the :func:`residuals`
    are i.i.d. unit-rate exponentials. Returns ``(D, p_value)`` from a one-sample KS
    test of the residuals against the ``Exp(1)`` CDF; a small p-value rejects the
    fitted model. Requires at least three events.
    """
    from .gof_tests import _ks_pvalue

    taus = residuals(events, mu, alpha, beta)
    m = len(taus)
    if m < 2:
        raise ValueError("need at least three events")
    s = sorted(taus)
    d = 0.0
    for i, x in enumerate(s):
        cdf = 1.0 - math.exp(-x)          # Exp(1) CDF
        d = max(d, abs((i + 1) / m - cdf), abs(cdf - i / m))
    return d, _ks_pvalue(d, math.sqrt(m))     # _ks_pvalue expects the effective size
