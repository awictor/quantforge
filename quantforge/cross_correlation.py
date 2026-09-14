"""FFT cross-correlation and lead-lag detection between two series.

The cross-correlation ``(x star y)[lag] = sum_n x[n+lag] y[n]`` measures how similar
two signals are when one is slid past the other. Its peak locates the *lag* at which
they line up best -- the standard way to estimate a time delay, a lead-lag relationship
between two price series, or the offset needed to align two recordings. Done directly
it costs ``O(n m)``; via the FFT (multiply one spectrum by the conjugate of the other)
it is ``O(N log N)``. Pure standard library on top of the radix-2 FFT.
"""

from .fft import fft, ifft


def _next_pow2(n):
    p = 1
    while p < n:
        p <<= 1
    return p


def cross_correlation(x, y, max_lag=None):
    """Cross-correlation of ``x`` and ``y`` at lags ``-max_lag..max_lag`` via the FFT.

    Returns ``(lags, values)`` where ``values[i]`` is the raw cross-correlation
    ``sum_n x[n + lags[i]] y[n]`` at ``lags[i]``. A positive lag slides ``x`` forward
    relative to ``y``; the peak lag is where the two series line up best. ``max_lag``
    defaults to ``len - 1``. Both inputs must be the same non-empty length.
    """
    n = len(x)
    if n == 0 or len(y) != n:
        raise ValueError("inputs must be non-empty and the same length")
    if max_lag is None:
        max_lag = n - 1
    if not (1 <= max_lag <= n - 1):
        raise ValueError("max_lag must be in [1, len-1]")
    # Zero-pad to >= 2n to make the circular FFT correlation act as a linear one.
    N = _next_pow2(2 * n)
    fx = fft(list(x) + [0.0] * (N - n))
    fy = fft(list(y) + [0.0] * (N - n))
    # X * conj(Y) -> correlation; bin k holds lag k (wrapping to negative lags at the top).
    corr = ifft([fx[k] * fy[k].conjugate() for k in range(N)])
    lags = list(range(-max_lag, max_lag + 1))
    values = []
    for lag in lags:
        idx = lag if lag >= 0 else N + lag
        values.append(corr[idx].real)
    return lags, values


def normalized_cross_correlation(x, y, max_lag=None):
    """Cross-correlation coefficient at each lag: dimensionless, in ``[-1, 1]``.

    Subtracts each series' mean and divides by ``sqrt(var_x * var_y) * n`` so the value
    is a correlation coefficient -- ``1`` at the lag of perfect alignment, ``0`` for
    uncorrelated series -- independent of the signals' amplitudes. Returns
    ``(lags, values)``.
    """
    n = len(x)
    if n == 0 or len(y) != n:
        raise ValueError("inputs must be non-empty and the same length")
    mx = sum(x) / n
    my = sum(y) / n
    xc = [v - mx for v in x]
    yc = [v - my for v in y]
    sx = sum(v * v for v in xc)
    sy = sum(v * v for v in yc)
    if sx <= 0.0 or sy <= 0.0:
        raise ValueError("zero-variance series")
    denom = (sx * sy) ** 0.5
    lags, raw = cross_correlation(xc, yc, max_lag)
    return lags, [v / denom for v in raw]


def lag_at_max_correlation(x, y, max_lag=None):
    """Lag (in samples) at which ``x`` and ``y`` are most positively correlated.

    The argmax of :func:`normalized_cross_correlation`. A positive result means ``x``
    *leads* ``y`` by that many samples (shift ``x`` forward to align them); a negative
    result means ``x`` lags ``y``. Also returns the peak correlation coefficient as
    ``(lag, coefficient)``.
    """
    lags, values = normalized_cross_correlation(x, y, max_lag)
    best = 0
    for i in range(1, len(values)):
        if values[i] > values[best]:
            best = i
    return lags[best], values[best]
