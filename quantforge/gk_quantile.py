"""Greenwald-Khanna quantile summary: deterministic epsilon-approximate quantiles.

The KLL sketch (:mod:`quantforge.kll`) achieves an optimal *randomized* rank-error bound. The
Greenwald-Khanna (2001) summary is the classic *deterministic* alternative: no coin flips, a
worst-case guarantee that any quantile query is answered within ``epsilon * n`` in rank, using
``O((1/epsilon) log(epsilon n))`` space. Each stored tuple keeps a value, ``g`` (the gap in rank
since the previous stored value) and ``delta`` (the uncertainty in that value's rank); an
insertion followed by periodic band-based compression keeps ``g + delta <= 2 epsilon n`` for
every tuple, which is exactly what bounds the query error.

Deterministic and reproducible with no seed. Pure standard library.
"""

import math


class GKQuantile:
    """Greenwald-Khanna epsilon-approximate quantile summary.

    ``epsilon`` is the maximum rank error as a fraction of ``n``. Feed values with
    :meth:`add`/:meth:`add_all`, query with :meth:`quantile` / :meth:`rank`. Any answer is within
    ``epsilon * n`` of the true rank.
    """

    def __init__(self, epsilon=0.01):
        if not 0 < epsilon < 1:
            raise ValueError("epsilon must be in (0, 1)")
        self.epsilon = float(epsilon)
        self._tuples = []     # list of [value, g, delta], sorted by value
        self._n = 0

    def add(self, x):
        """Insert one value."""
        x = float(x)
        # find insertion index (first tuple with value >= x)
        lo, hi = 0, len(self._tuples)
        while lo < hi:
            mid = (lo + hi) // 2
            if self._tuples[mid][0] < x:
                lo = mid + 1
            else:
                hi = mid
        i = lo
        if i == 0 or i == len(self._tuples):
            delta = 0                       # new min or max is exact
        else:
            delta = int(math.floor(2.0 * self.epsilon * self._n))
        self._tuples.insert(i, [x, 1, delta])
        self._n += 1
        if self._n % max(1, int(1.0 / (2.0 * self.epsilon))) == 0:
            self._compress()
        return self

    def add_all(self, xs):
        for x in xs:
            self.add(x)
        return self

    def _compress(self):
        cap = int(2.0 * self.epsilon * self._n)
        i = len(self._tuples) - 2
        while i >= 1:
            # merge tuple i into i+1 if the combined band allows it
            if (self._tuples[i][1] + self._tuples[i + 1][1] + self._tuples[i + 1][2]) <= cap:
                self._tuples[i + 1][1] += self._tuples[i][1]
                del self._tuples[i]
            i -= 1

    @property
    def n(self):
        return self._n

    def rank(self, x):
        """Estimated number of elements ``<= x`` (mid-band estimate)."""
        r = 0
        for v, g, delta in self._tuples:
            if v <= x:
                r += g
            else:
                break
        return r

    def quantile(self, q):
        """Approximate value at quantile ``q`` in ``[0, 1]`` within ``epsilon * n`` rank error."""
        if self._n == 0:
            raise ValueError("empty summary")
        target = q * self._n
        allowed = self.epsilon * self._n
        rmin = 0
        best = self._tuples[-1][0]
        for v, g, delta in self._tuples:
            rmin += g                       # rmin_i = cumulative g through tuple i
            rmax = rmin + delta             # rmax_i = rmin_i + delta_i
            # Greenwald-Khanna query: tuple i answers if its rank band brackets target
            # within the error on both sides.
            if target - rmin <= allowed and rmax - target <= allowed:
                return v
            if rmin > target + allowed:
                break
        return best
