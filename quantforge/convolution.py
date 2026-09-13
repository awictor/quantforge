"""FFT-based linear convolution and autocorrelation.

Linear convolution ``(a * b)[k] = sum_i a[i] b[k-i]`` is the coefficient product of
two polynomials and the workhorse of filtering. Done directly it costs ``O(n m)``;
via the FFT it is ``O(N log N)`` after zero-padding both inputs to a common
power-of-two length ``N >= len(a)+len(b)-1`` and multiplying their spectra. The same
trick gives the full autocorrelation in ``O(n log n)``. Pure standard library on top
of the radix-2 FFT.
"""

from .fft import fft, ifft


def _next_pow2(n):
    p = 1
    while p < n:
        p <<= 1
    return p


def convolve(a, b):
    """Full linear convolution of ``a`` and ``b`` via the FFT.

    Returns a real list of length ``len(a) + len(b) - 1`` -- equivalently the
    coefficients of the product of the two polynomials with coefficients ``a`` and
    ``b``. ``O(N log N)`` where ``N`` is the padded power-of-two length.
    """
    if not a or not b:
        raise ValueError("inputs must be non-empty")
    out_len = len(a) + len(b) - 1
    n = _next_pow2(out_len)
    fa = fft(list(a) + [0.0] * (n - len(a)))
    fb = fft(list(b) + [0.0] * (n - len(b)))
    prod = [fa[i] * fb[i] for i in range(n)]
    conv = ifft(prod)
    return [conv[i].real for i in range(out_len)]


def fft_autocorrelation(x, max_lag=None):
    """Autocorrelation of ``x`` at lags ``0..max_lag`` via the FFT.

    Uses the Wiener-Khinchin route (spectrum times its conjugate) on the
    mean-subtracted series, then normalizes by lag-0 so ``acf[0] = 1``. ``max_lag``
    defaults to ``len(x) - 1``. Matches a direct autocovariance sum but in
    ``O(n log n)``.
    """
    n = len(x)
    if n < 2:
        raise ValueError("need at least two points")
    if max_lag is None:
        max_lag = n - 1
    mean = sum(x) / n
    xc = [v - mean for v in x]
    # Zero-pad to >= 2n to avoid circular wrap-around.
    N = _next_pow2(2 * n)
    fx = fft(xc + [0.0] * (N - n))
    power = [fx[i] * fx[i].conjugate() for i in range(N)]
    acov = ifft(power)
    c0 = acov[0].real
    if c0 <= 0.0:
        raise ValueError("zero-variance series")
    return [acov[k].real / c0 for k in range(max_lag + 1)]
