"""Radix-2 Cooley-Tukey fast Fourier transform (public API).

``fft`` computes the discrete Fourier transform of a complex (or real) sequence in
``O(n log n)`` instead of the ``O(n^2)`` direct sum, and ``ifft`` inverts it. The
length must be a power of two. These wrap the same radix-2 core used by the
Carr-Madan option-pricing strip, exposed for general spectral work. Pure standard
library.
"""

import cmath
import math


def _radix2(x, inverse):
    n = len(x)
    if n == 0 or (n & (n - 1)) != 0:
        raise ValueError("length must be a positive power of two")
    a = [complex(v) for v in x]
    # Bit-reversal permutation.
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    sign = 1.0 if inverse else -1.0
    while length <= n:
        wlen = cmath.exp(1j * sign * 2.0 * math.pi / length)
        half = length >> 1
        for start in range(0, n, length):
            w = 1.0 + 0j
            for k in range(half):
                u = a[start + k]
                v = a[start + k + half] * w
                a[start + k] = u + v
                a[start + k + half] = u - v
                w *= wlen
        length <<= 1
    if inverse:
        a = [v / n for v in a]
    return a


def fft(x):
    """Forward FFT of a sequence whose length is a power of two.

    Returns the complex DFT ``X_k = sum_n x_n exp(-2 pi i k n / N)``. Accepts real or
    complex input; raises unless the length is a positive power of two.
    """
    return _radix2(x, inverse=False)


def ifft(x):
    """Inverse FFT: recovers the sequence from its DFT (with the ``1/N`` scaling).

    ``ifft(fft(x)) == x`` up to floating error. Length must be a power of two.
    """
    return _radix2(x, inverse=True)
