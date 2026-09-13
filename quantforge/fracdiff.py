"""Fractional differencing for long-memory (ARFIMA) time series.

Integer differencing (``d = 1``) removes a unit root but wipes out the memory in a
series. Fractional differencing applies the operator ``(1 - L)^d`` for real ``d`` via
its binomial expansion,

    (1 - L)^d = sum_{k=0}^inf w_k L^k,   w_0 = 1,  w_k = w_{k-1} * -(d - k + 1) / k,

so a small ``d`` can make a series stationary while keeping most of its
autocorrelation (Hosking 1981). The weights decay slowly, so two truncations are
offered: a full expansion up to the series length, and Lopez de Prado's
fixed-width window that drops weights below a tolerance. Pure standard library.
"""


def fracdiff_weights(d, n):
    """Binomial weights ``w_0..w_{n-1}`` of the operator ``(1 - L)^d``.

    ``w_0 = 1`` and ``w_k = w_{k-1} * -(d - k + 1) / k``. For integer ``d`` the
    weights vanish beyond ``k = d`` (e.g. ``d = 1`` gives ``[1, -1]`` then zeros).
    """
    if n < 1:
        raise ValueError("n must be at least 1")
    w = [1.0]
    for k in range(1, n):
        w.append(w[-1] * -(d - k + 1) / k)
    return w


def fractional_difference(series, d):
    """Fractionally difference ``series`` by the full expansion of ``(1 - L)^d``.

    Returns a list the same length as ``series``; entry ``t`` is
    ``sum_{k=0}^{t} w_k * series[t - k]`` (a growing backward window, so early
    entries use fewer weights). ``d = 0`` returns the series unchanged and ``d = 1``
    returns the first difference (with the first entry equal to ``series[0]``).
    """
    n = len(series)
    if n == 0:
        raise ValueError("series must be non-empty")
    w = fracdiff_weights(d, n)
    out = []
    for t in range(n):
        out.append(sum(w[k] * series[t - k] for k in range(t + 1)))
    return out


def fixed_width_fracdiff(series, d, threshold=1e-5):
    """Fixed-width fractional differencing (Lopez de Prado).

    Truncates the weight window at the first ``|w_k| < threshold`` and applies that
    fixed-length filter, so every output point uses the same weights. Returns the
    valid portion (length ``len(series) - width + 1``), where ``width`` is the number
    of retained weights. Preserves stationarity with a constant memory window.
    """
    n = len(series)
    if n == 0:
        raise ValueError("series must be non-empty")
    # Grow weights until they fall below the threshold.
    w = [1.0]
    k = 1
    while k < n:
        wk = w[-1] * -(d - k + 1) / k
        if abs(wk) < threshold:
            break
        w.append(wk)
        k += 1
    width = len(w)
    if width > n:
        raise ValueError("threshold too small for series length")
    out = []
    for t in range(width - 1, n):
        out.append(sum(w[j] * series[t - j] for j in range(width)))
    return out
