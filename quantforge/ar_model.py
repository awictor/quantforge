"""Autoregressive model fitting by Yule-Walker (Durbin-Levinson).

An AR(p) process ``x_t = c + sum_{i=1}^p phi_i x_{t-i} + eps_t`` is fit by solving
the Yule-Walker equations, which relate the AR coefficients to the
autocovariances. The Durbin-Levinson recursion solves them in O(p^2) and also
yields the innovation (white-noise) variance at each order. The intercept is
recovered from the process mean: ``c = mean * (1 - sum phi)``. Pure standard
library.
"""


def _autocov(x, k, mean):
    n = len(x)
    return sum((x[t] - mean) * (x[t - k] - mean) for t in range(k, n)) / n


def fit_ar_yule_walker(x, order):
    """Fit an AR(``order``) model by the Yule-Walker / Durbin-Levinson method.

    Parameters
    ----------
    x : sequence of float
        The series.
    order : int
        AR order ``p`` (>= 1).

    Returns
    -------
    dict
        ``{"coefficients": [phi_1, ..., phi_p], "intercept": c,
        "noise_variance": sigma2, "mean": mu}``. The coefficients solve the
        Yule-Walker equations; ``noise_variance`` is the innovation variance from
        the final Levinson step.
    """
    n = len(x)
    if order < 1:
        raise ValueError("order must be >= 1")
    if n <= order + 1:
        raise ValueError("series too short for the requested order")
    mean = sum(x) / n
    g = [_autocov(x, k, mean) for k in range(order + 1)]
    if g[0] == 0.0:
        raise ValueError("series has zero variance")
    r = [g[k] / g[0] for k in range(order + 1)]   # autocorrelations

    # Durbin-Levinson recursion.
    phi = [0.0] * (order + 1)
    phi[1] = r[1]
    v = 1.0 - r[1] * r[1]                          # normalized innovation variance
    for k in range(2, order + 1):
        acc = r[k] - sum(phi[j] * r[k - j] for j in range(1, k))
        phi_kk = acc / v if v != 0.0 else 0.0
        new = phi[:]
        for j in range(1, k):
            new[j] = phi[j] - phi_kk * phi[k - j]
        new[k] = phi_kk
        phi = new
        v *= (1.0 - phi_kk * phi_kk)

    coeffs = phi[1:order + 1]
    intercept = mean * (1.0 - sum(coeffs))
    noise_var = v * g[0]
    return {
        "coefficients": coeffs,
        "intercept": intercept,
        "noise_variance": noise_var,
        "mean": mean,
    }


def ar_forecast(model, history, steps=1):
    """Forecast ``steps`` ahead from a fitted AR model and recent ``history``.

    ``history`` must hold at least ``len(coefficients)`` most-recent observations
    (oldest first). Iterates the deterministic AR recursion (innovations set to
    their zero mean), appending each forecast to drive the next.
    """
    coeffs = model["coefficients"]
    c = model["intercept"]
    p = len(coeffs)
    if len(history) < p:
        raise ValueError("history shorter than the AR order")
    if steps < 1:
        raise ValueError("steps must be >= 1")
    series = list(history)
    out = []
    for _ in range(steps):
        # phi_1 multiplies the most recent point.
        pred = c + sum(coeffs[i] * series[-1 - i] for i in range(p))
        series.append(pred)
        out.append(pred)
    return out
