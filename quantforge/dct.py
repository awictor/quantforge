"""Discrete cosine transform (DCT-II) and its inverse (DCT-III).

The DCT expresses a real signal as a sum of cosines at different frequencies -- the
real-valued cousin of the FFT and the backbone of JPEG/MP3 compression, because it packs
most of a smooth signal's energy into the first few coefficients. This implements the
orthonormal DCT-II (forward) and DCT-III (inverse) so that ``idct(dct(x)) == x`` and the
transform preserves energy (Parseval). Direct O(n^2) evaluation. Pure standard library.
"""

import math


def dct(x):
    """Orthonormal DCT-II of a real sequence.

    ``X[k] = s(k) sum_n x[n] cos(pi (2n+1) k / (2N))`` with the orthonormal scaling
    ``s(0) = sqrt(1/N)``, ``s(k>0) = sqrt(2/N)``. Energy-preserving; pair with
    :func:`idct`.
    """
    n = len(x)
    if n == 0:
        raise ValueError("need at least one sample")
    out = []
    for k in range(n):
        s = math.sqrt(1.0 / n) if k == 0 else math.sqrt(2.0 / n)
        total = sum(x[i] * math.cos(math.pi * (2 * i + 1) * k / (2 * n))
                    for i in range(n))
        out.append(s * total)
    return out


def idct(X):
    """Orthonormal inverse DCT (DCT-III) -- exact inverse of :func:`dct`.

    ``x[n] = sum_k s(k) X[k] cos(pi (2n+1) k / (2N))`` with the same orthonormal scaling.
    Recovers the original sequence to machine precision.
    """
    n = len(X)
    if n == 0:
        raise ValueError("need at least one coefficient")
    out = []
    for i in range(n):
        total = 0.0
        for k in range(n):
            s = math.sqrt(1.0 / n) if k == 0 else math.sqrt(2.0 / n)
            total += s * X[k] * math.cos(math.pi * (2 * i + 1) * k / (2 * n))
        out.append(total)
    return out
