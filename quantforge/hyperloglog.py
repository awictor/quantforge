"""HyperLogLog: streaming distinct-count (cardinality) estimation.

Counting *distinct* items exactly needs memory proportional to the number of distinct
items. HyperLogLog (Flajolet et al. 2007) estimates the count from a fixed number of
small registers -- a few kilobytes handles billions of distinct items with a typical
error of about ``1.04 / sqrt(m)``. Each item is hashed; the first ``p`` bits pick a
register, and each register keeps the largest number of leading zeros seen in the rest
of the hash. The harmonic mean of ``2^register`` across registers, times a bias
constant, estimates the cardinality (with small- and large-range corrections). Uses
``hashlib`` for a stable hash. Pure standard library.
"""

import hashlib
import math


class HyperLogLog:
    """HyperLogLog distinct-count estimator with ``2^p`` registers.

    ``p`` in ``[4, 16]`` trades memory for accuracy: ``m = 2^p`` registers give a
    relative error near ``1.04 / sqrt(m)`` (``p = 14`` -> ~0.8%). Add items with
    :meth:`add` (any hashable stringifiable value); read the estimate from
    :meth:`count`. Two sketches with the same ``p`` merge via :meth:`merge` (union
    cardinality), the property that makes it distributable.
    """

    def __init__(self, p=14):
        if not (4 <= p <= 16):
            raise ValueError("p must be in [4, 16]")
        self.p = p
        self.m = 1 << p
        self.registers = [0] * self.m
        self.alpha = self._alpha(self.m)

    @staticmethod
    def _alpha(m):
        if m == 16:
            return 0.673
        if m == 32:
            return 0.697
        if m == 64:
            return 0.709
        return 0.7213 / (1.0 + 1.079 / m)

    def _hash(self, value):
        h = hashlib.sha1(str(value).encode("utf-8")).digest()
        # Take 64 bits.
        return int.from_bytes(h[:8], "big")

    def add(self, value):
        """Add an item to the sketch."""
        x = self._hash(value)
        idx = x >> (64 - self.p)                 # top p bits -> register index
        rest = (x << self.p) & ((1 << 64) - 1)   # remaining bits, left-aligned in 64
        # Position (1-based) of the leftmost set bit in the remaining (64-p) bits.
        rank = 1
        mask = 1 << 63
        while rank <= 64 - self.p and not (rest & mask):
            rank += 1
            mask >>= 1
        if rank > self.registers[idx]:
            self.registers[idx] = rank
        return self

    def count(self):
        """Estimate the number of distinct items added."""
        m = self.m
        raw = self.alpha * m * m / sum(2.0 ** (-r) for r in self.registers)
        if raw <= 2.5 * m:
            # Small-range correction: linear counting on empty registers.
            zeros = self.registers.count(0)
            if zeros != 0:
                return m * math.log(m / zeros)
        return raw

    def merge(self, other):
        """Merge ``other`` into this sketch (element-wise max); returns self."""
        if other.p != self.p:
            raise ValueError("cannot merge sketches with different p")
        self.registers = [max(a, b) for a, b in zip(self.registers, other.registers)]
        return self
