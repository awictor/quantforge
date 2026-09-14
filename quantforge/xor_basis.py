"""XOR linear basis over GF(2): the vector space spanned by a set of integers under XOR.

Treating each non-negative integer as a bit vector over GF(2), a set of integers spans a
linear space under XOR. Gaussian elimination keeps one basis vector per leading bit, which
answers the standard questions in ``O(bits)`` each: the maximum (or minimum) XOR of any
subset, whether a value is representable, how many distinct XOR values exist (``2^rank``),
and the ``k``-th smallest reachable value. Inserting ``n`` numbers of ``b`` bits costs
``O(n * b)``. Pure standard library.
"""


class XorBasis:
    """A GF(2) linear basis maintained by leading bit.

    ``insert(x)`` adds an integer to the span (reducing it against the current basis);
    ``max_xor``/``min_xor`` return the extremal subset-XOR, ``can_represent`` tests
    membership, ``rank`` is the basis size, ``count_distinct`` is ``2^rank``, and
    ``kth_smallest`` indexes the sorted reachable values. Non-negative integers only.
    """

    __slots__ = ("_basis", "_size")

    def __init__(self, values=None):
        self._basis = {}          # leading-bit position -> basis vector
        self._size = 0
        if values:
            for v in values:
                self.insert(v)

    def insert(self, x):
        """Add ``x`` to the span. Returns ``True`` if it enlarged the basis (was independent)."""
        if x < 0:
            raise ValueError("XOR basis handles non-negative integers only")
        cur = x
        while cur:
            hb = cur.bit_length() - 1
            if hb not in self._basis:
                self._basis[hb] = cur
                self._size += 1
                return True
            cur ^= self._basis[hb]
        return False              # reduced to zero: already representable

    def __len__(self):
        return self._size

    def rank(self):
        """Number of independent basis vectors (dimension of the span)."""
        return self._size

    def can_represent(self, x):
        """True if ``x`` is the XOR of some subset of the inserted values."""
        if x < 0:
            return False
        cur = x
        while cur:
            hb = cur.bit_length() - 1
            if hb not in self._basis:
                return False
            cur ^= self._basis[hb]
        return True

    def max_xor(self, start=0):
        """Maximum value of ``start`` XOR any subset of the span (default ``start=0``)."""
        res = start
        for hb in sorted(self._basis, reverse=True):
            if res ^ self._basis[hb] > res:
                res ^= self._basis[hb]
        return res

    def min_xor(self, start=0):
        """Minimum value of ``start`` XOR any subset of the span."""
        res = start
        for hb in sorted(self._basis, reverse=True):
            if res ^ self._basis[hb] < res:
                res ^= self._basis[hb]
        return res

    def count_distinct(self):
        """Number of distinct XOR values reachable, i.e. ``2^rank`` (includes 0)."""
        return 1 << self._size

    def kth_smallest(self, k):
        """The ``k``-th smallest distinct reachable XOR value (``0``-indexed; ``0`` -> 0).

        Reduces the basis to a canonical form where each vector owns a distinct leading bit
        with zeros elsewhere among the pivots, so the ``k``-th value is the XOR of the
        pivots picked out by the bits of ``k``.
        """
        if not 0 <= k < self.count_distinct():
            raise IndexError("k out of range [0, 2^rank)")
        # canonicalize: reduce each pivot so no other pivot shares its lower bits
        pivots = sorted(self._basis)          # ascending leading-bit order
        reduced = [self._basis[hb] for hb in pivots]
        for i in range(len(reduced)):
            for j in range(i):
                if reduced[i] ^ reduced[j] < reduced[i]:
                    reduced[i] ^= reduced[j]
        result = 0
        for i in range(len(reduced)):
            if (k >> i) & 1:
                result ^= reduced[i]
        return result

    def merge(self, other):
        """Fold another basis's span into this one."""
        if not isinstance(other, XorBasis):
            raise TypeError("can only merge with another XorBasis")
        for vec in other._basis.values():
            self.insert(vec)
        return self
