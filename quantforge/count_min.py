"""Count-Min sketch and Bloom filter: sublinear frequency and membership.

Two probabilistic structures for streams too large to store exactly:

  * ``CountMinSketch`` -- estimates the frequency of any item from a fixed ``depth x
    width`` table of counters. Each item is hashed into one counter per row and all are
    incremented; the estimate is the *minimum* across rows, which can only
    *over*-count (collisions add, never subtract), with the error bounded by the table
    width.
  * ``BloomFilter`` -- tests set membership from a bit array. It never reports a false
    negative (a member always tests positive) but allows a tunable false-positive rate.

Both use ``hashlib`` for stable, seedable hashing. Pure standard library.
"""

import hashlib
import math


def _hashes(value, k, mod):
    """``k`` independent hash values in ``[0, mod)`` from one SHA-1 digest (double hashing)."""
    data = str(value).encode("utf-8")
    h1 = int.from_bytes(hashlib.sha1(data + b"\x00").digest()[:8], "big")
    h2 = int.from_bytes(hashlib.sha1(data + b"\x01").digest()[:8], "big")
    return [(h1 + i * h2) % mod for i in range(k)]


class CountMinSketch:
    """Count-Min sketch for streaming frequency estimation.

    ``width`` counters per row, ``depth`` rows. The estimate never under-counts; the
    over-count is at most about ``total_added / width`` with probability
    ``1 - (1/2)^depth``. Add occurrences with :meth:`add` (optionally a count), query
    with :meth:`estimate`.
    """

    def __init__(self, width=1024, depth=5):
        if width < 1 or depth < 1:
            raise ValueError("width and depth must be >= 1")
        self.width = width
        self.depth = depth
        self.table = [[0] * width for _ in range(depth)]
        self.total = 0

    def add(self, value, count=1):
        """Add ``count`` occurrences of ``value``."""
        hs = _hashes(value, self.depth, self.width)
        for r in range(self.depth):
            self.table[r][hs[r]] += count
        self.total += count
        return self

    def estimate(self, value):
        """Estimated frequency of ``value`` (an upper bound on the true count)."""
        hs = _hashes(value, self.depth, self.width)
        return min(self.table[r][hs[r]] for r in range(self.depth))


class BloomFilter:
    """Bloom filter for approximate set membership (no false negatives).

    ``capacity`` expected items and target ``error_rate`` size the bit array and hash
    count optimally. :meth:`add` inserts, :meth:`contains` tests: a member always
    returns True; a non-member returns True only with probability ~``error_rate``.
    """

    def __init__(self, capacity=1000, error_rate=0.01):
        if capacity < 1:
            raise ValueError("capacity must be >= 1")
        if not (0.0 < error_rate < 1.0):
            raise ValueError("error_rate must be in (0, 1)")
        # Optimal bit count m and hash count k.
        self.m = max(1, int(-capacity * math.log(error_rate) / (math.log(2) ** 2)))
        self.k = max(1, round(self.m / capacity * math.log(2)))
        self.bits = bytearray((self.m + 7) // 8)
        self.n = 0

    def _set(self, i):
        self.bits[i >> 3] |= (1 << (i & 7))

    def _get(self, i):
        return (self.bits[i >> 3] >> (i & 7)) & 1

    def add(self, value):
        """Insert ``value`` into the filter."""
        for i in _hashes(value, self.k, self.m):
            self._set(i)
        self.n += 1
        return self

    def contains(self, value):
        """Test membership: True if possibly present, False if definitely absent."""
        return all(self._get(i) for i in _hashes(value, self.k, self.m))

    def __contains__(self, value):
        return self.contains(value)
