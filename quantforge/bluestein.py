"""Bluestein's algorithm: the DFT of an arbitrary-length sequence.

The radix-2 FFT in :mod:`quantforge.fft` needs a power-of-two length. Bluestein's algorithm
(the chirp-Z transform) computes the DFT of *any* length ``n`` in ``O(n log n)`` by
rewriting each output as a convolution: ``k n = (k^2 + n^2 - (k-n)^2) / 2`` turns the DFT
into a multiplication by a chirp, a convolution with another chirp, and a final chirp
multiply. The convolution is padded up to a power of two and done with the existing radix-2
FFT. Pure standard library (uses `cmath`).
"""

import cmath

from .fft import fft as _fft, ifft as _ifft


def _next_pow2(m):
    p = 1
    while p < m:
        p <<= 1
    return p


def _bluestein(x, inverse):
    n = len(x)
    if n == 0:
        return []
    if n == 1:
        return [complex(x[0])]
    sign = 1.0 if inverse else -1.0
    # chirp a_k = x_k * exp(sign * pi i k^2 / n)
    w = [cmath.exp(sign * 1j * cmath.pi * (k * k % (2 * n)) / n) for k in range(n)]
    a = [complex(x[k]) * w[k] for k in range(n)]
    # b_k = conj(w_k), extended symmetrically for the linear convolution
    m = _next_pow2(2 * n - 1)
    a_pad = a + [0j] * (m - n)
    b_pad = [0j] * m
    for k in range(n):
        b_pad[k] = w[k].conjugate()
        if k:
            b_pad[m - k] = w[k].conjugate()
    # circular convolution via FFT
    fa = _fft(a_pad)
    fb = _fft(b_pad)
    fc = [fa[i] * fb[i] for i in range(m)]
    conv = _ifft(fc)
    out = [conv[k] * w[k] for k in range(n)]
    if inverse:
        out = [v / n for v in out]
    return out


def dft(x):
    """Forward DFT of a sequence of *any* length ``n`` (Bluestein), returning complex output.

    ``X_k = sum_n x_n exp(-2 pi i k n / N)``. Matches :func:`quantforge.fft.fft` on
    power-of-two lengths and a direct DFT everywhere, in ``O(n log n)``.
    """
    return _bluestein(list(x), inverse=False)


def idft(x):
    """Inverse DFT of any length (Bluestein), with the ``1/N`` scaling.

    ``idft(dft(x)) == x`` up to floating error.
    """
    return _bluestein(list(x), inverse=True)
