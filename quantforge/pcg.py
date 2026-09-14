"""Statistically strong pseudo-random generators: PCG32 and xorshift128+.

The library's other stochastic code uses a plain LCG for reproducibility; its low bits are
weak. PCG32 (O'Neill's permuted congruential generator) and xorshift128+ pass the
empirical tests an LCG fails, while staying tiny and fast. Both produce reproducible
streams from a seed and expose uniform ``[0, 1)`` floats. Pure standard library.
"""

_MASK64 = (1 << 64) - 1
_MASK32 = (1 << 32) - 1


class PCG32:
    """PCG-XSH-RR 32-bit generator (O'Neill). Strong, small, reproducible.

    ``next_uint32()`` yields a 32-bit output; ``random()`` a float in ``[0, 1)``;
    ``randint(lo, hi)`` an integer in ``[lo, hi]`` (inclusive) without modulo bias.
    Seeded by ``seed`` and an optional stream ``seq``.
    """

    _MULT = 6364136223846793005

    def __init__(self, seed=1234567, seq=54):
        self.inc = ((seq << 1) | 1) & _MASK64
        self.state = 0
        self.next_uint32()
        self.state = (self.state + (seed & _MASK64)) & _MASK64
        self.next_uint32()

    def next_uint32(self):
        """Advance the state and return the next 32-bit output."""
        old = self.state
        self.state = (old * self._MULT + self.inc) & _MASK64
        xorshifted = (((old >> 18) ^ old) >> 27) & _MASK32
        rot = old >> 59
        return ((xorshifted >> rot) | (xorshifted << ((-rot) & 31))) & _MASK32

    def random(self):
        """A uniform float in ``[0, 1)`` (24-bit-plus resolution from one 32-bit word)."""
        return self.next_uint32() / 4294967296.0

    def randint(self, lo, hi):
        """A uniform integer in ``[lo, hi]`` inclusive, rejection-sampled (no modulo bias)."""
        if hi < lo:
            raise ValueError("require lo <= hi")
        span = hi - lo + 1
        # Rejection threshold to remove the bias of the top partial bucket.
        threshold = (0x100000000 % span)
        while True:
            r = self.next_uint32()
            if r >= threshold:
                return lo + (r % span)


class Xorshift128Plus:
    """xorshift128+ generator (Vigna). 64-bit output, long period, reproducible."""

    def __init__(self, seed=1234567):
        # Seed the two 64-bit words via a splitmix64 warm-up (avoids all-zero state).
        s = seed & _MASK64
        self.s0 = self._splitmix(s)
        self.s1 = self._splitmix(self.s0)
        if self.s0 == 0 and self.s1 == 0:
            self.s0 = 1

    @staticmethod
    def _splitmix(x):
        x = (x + 0x9E3779B97F4A7C15) & _MASK64
        z = x
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _MASK64
        return (z ^ (z >> 31)) & _MASK64

    def next_uint64(self):
        """Advance and return the next 64-bit output."""
        x = self.s0
        y = self.s1
        self.s0 = y
        x ^= (x << 23) & _MASK64
        self.s1 = (x ^ y ^ (x >> 17) ^ (y >> 26)) & _MASK64
        return (self.s1 + y) & _MASK64

    def random(self):
        """A uniform float in ``[0, 1)`` from the top 53 bits of a 64-bit output."""
        return (self.next_uint64() >> 11) / 9007199254740992.0
