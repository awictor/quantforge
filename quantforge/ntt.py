"""Number-theoretic transform (NTT): exact integer convolution modulo a prime.

The NTT is the FFT done in a finite field instead of the complex numbers: it replaces the
complex root of unity ``e^{-2 pi i / n}`` with a primitive ``n``-th root of unity modulo a
prime, so it carries *no floating-point error*. With the standard prime
``998244353 = 119 * 2^23 + 1`` it supports transform lengths up to ``2^23``. Its headline
use is exact convolution -- multiplying two integer polynomials (or big integers) with all
coefficients computed mod the prime -- in ``O(n log n)``. Pure standard library.
"""

# 998244353 = 119 * 2^23 + 1; 3 is a primitive root.
NTT_PRIME = 998244353
_PRIMITIVE_ROOT = 3


def _next_pow2(m):
    p = 1
    while p < m:
        p <<= 1
    return p


def _transform(a, invert, mod, root):
    n = len(a)
    a = a[:]
    # bit-reversal permutation
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    while length <= n:
        # w_len is a primitive length-th root of unity mod `mod`
        if invert:
            wlen = pow(root, mod - 1 - (mod - 1) // length, mod)
        else:
            wlen = pow(root, (mod - 1) // length, mod)
        half = length >> 1
        for start in range(0, n, length):
            w = 1
            for k in range(half):
                u = a[start + k]
                v = a[start + k + half] * w % mod
                a[start + k] = (u + v) % mod
                a[start + k + half] = (u - v) % mod
                w = w * wlen % mod
        length <<= 1
    if invert:
        inv_n = pow(n, mod - 2, mod)
        a = [x * inv_n % mod for x in a]
    return a


def ntt(a, mod=NTT_PRIME):
    """Forward number-theoretic transform of ``a`` (length must be a power of two)."""
    n = len(a)
    if n & (n - 1) != 0:
        raise ValueError("length must be a power of two")
    return _transform([x % mod for x in a], invert=False, mod=mod, root=_PRIMITIVE_ROOT)


def intt(a, mod=NTT_PRIME):
    """Inverse number-theoretic transform (with the ``1/n`` field scaling)."""
    n = len(a)
    if n & (n - 1) != 0:
        raise ValueError("length must be a power of two")
    return _transform([x % mod for x in a], invert=True, mod=mod, root=_PRIMITIVE_ROOT)


def convolve_mod(a, b, mod=NTT_PRIME):
    """Exact convolution (polynomial product) of integer sequences ``a`` and ``b`` mod ``mod``.

    Returns the ``len(a) + len(b) - 1`` coefficients of the product, each reduced mod
    ``mod`` -- no floating-point error. Empty inputs give an empty result.
    """
    if not a or not b:
        return []
    result_len = len(a) + len(b) - 1
    size = _next_pow2(result_len)
    fa = _transform([x % mod for x in a] + [0] * (size - len(a)), False, mod, _PRIMITIVE_ROOT)
    fb = _transform([x % mod for x in b] + [0] * (size - len(b)), False, mod, _PRIMITIVE_ROOT)
    fc = [fa[i] * fb[i] % mod for i in range(size)]
    conv = _transform(fc, True, mod, _PRIMITIVE_ROOT)
    return conv[:result_len]
