"""Real cepstrum and quefrency analysis (echo and pitch detection).

The cepstrum is the inverse transform of the log-magnitude spectrum. Where the spectrum
turns periodicity in time into peaks in frequency, the cepstrum turns *periodicity in
the spectrum* -- the evenly spaced harmonics of a pitched sound, or the ripple a single
echo adds -- into a peak at a single "quefrency" (a time, in samples). That makes it the
classic tool for estimating the fundamental period of a voiced signal and for detecting
the delay of an echo. Pure standard library on top of the radix-2 FFT.
"""

import cmath
import math

from .fft import fft, ifft


def _next_pow2(n):
    p = 1
    while p < n:
        p <<= 1
    return p


def real_cepstrum(x):
    """Real cepstrum of ``x``: ``ifft(log|fft(x)|)``, real part.

    The inverse FFT of the log-magnitude spectrum. Returned length equals the padded
    power-of-two length used internally (``>= len(x)``). The independent variable is
    *quefrency* -- an index in samples, read like a lag. A tiny floor is added inside
    the log to keep spectral nulls finite.
    """
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    N = _next_pow2(n)
    X = fft(list(x) + [0.0] * (N - n))
    floor = 1e-300
    log_mag = [complex(math.log(abs(v) + floor), 0.0) for v in X]
    c = ifft(log_mag)
    return [v.real for v in c]


def power_cepstrum(x):
    """Power cepstrum of ``x``: ``|ifft(log|fft(x)|^2)|^2``.

    Squares the log-magnitude spectrum before inverting, then takes the squared
    magnitude -- emphasizing quefrency peaks. Same length convention as
    :func:`real_cepstrum`.
    """
    n = len(x)
    if n == 0:
        raise ValueError("input must be non-empty")
    N = _next_pow2(n)
    X = fft(list(x) + [0.0] * (N - n))
    floor = 1e-300
    log_pow = [complex(2.0 * math.log(abs(v) + floor), 0.0) for v in X]
    c = ifft(log_pow)
    return [abs(v) ** 2 for v in c]


def fundamental_quefrency(x, min_quefrency=1, max_quefrency=None):
    """Quefrency (in samples) of the largest real-cepstrum peak in a search band.

    For a voiced/pitched or echoed signal this is the fundamental period (or echo
    delay) in samples: divide the sample rate by it to get the pitch in Hz. The search
    ignores the low-quefrency region below ``min_quefrency`` (the spectral envelope) and
    is capped at ``max_quefrency`` (default: half the cepstrum length, the useful range
    of a real cepstrum). Returns ``(quefrency, peak_value)``.
    """
    c = real_cepstrum(x)
    n = len(c)
    if max_quefrency is None:
        max_quefrency = n // 2
    lo = max(1, min_quefrency)
    hi = min(max_quefrency, n // 2)
    if lo > hi:
        raise ValueError("empty quefrency search band")
    best = lo
    for q in range(lo + 1, hi + 1):
        if c[q] > c[best]:
            best = q
    return best, c[best]
