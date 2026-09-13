"""Streaming quantile estimation (P-square) and reservoir sampling.

Both process a stream in one pass with bounded memory -- no storing or sorting the data:

  * ``P2Quantile`` -- Jain & Chlamtac's P-square algorithm tracks a single quantile
    (e.g. the median or the 95th percentile) using just five markers whose positions
    adjust as data arrives, updating them by a piecewise-parabolic rule. Memory is
    O(1) regardless of stream length.
  * ``reservoir_sample`` -- Vitter's algorithm R draws a uniform sample of ``k`` items
    from a stream of unknown length in one pass, each item equally likely to be kept.

Pure standard library; the reservoir sampler uses a seeded LCG for reproducibility.
"""


class P2Quantile:
    """P-square streaming estimator of a single quantile ``p`` in ``[0, 1]``.

    Feed values with :meth:`update`; read the current estimate from :meth:`value`.
    Uses five markers and O(1) memory; accuracy improves as the stream grows. Best for
    smooth, stationary streams -- it approximates the true quantile without storing data.
    """

    def __init__(self, p):
        if not (0.0 < p < 1.0):
            raise ValueError("p must be in (0, 1)")
        self.p = p
        self.n = 0
        self.q = []        # marker heights
        self.npos = []     # marker positions (1-based)
        self.nprime = []   # desired positions
        self.dn = [0.0, p / 2.0, p, (1.0 + p) / 2.0, 1.0]

    def update(self, x):
        """Incorporate one observation ``x``."""
        if self.n < 5:
            self.q.append(x)
            self.n += 1
            if self.n == 5:
                self.q.sort()
                self.npos = [1, 2, 3, 4, 5]
                self.nprime = [1.0,
                               1.0 + 2.0 * self.p,
                               1.0 + 4.0 * self.p,
                               3.0 + 2.0 * self.p,
                               5.0]
            return self

        # Find cell k.
        if x < self.q[0]:
            self.q[0] = x
            k = 0
        elif x >= self.q[4]:
            self.q[4] = x
            k = 3
        else:
            k = 0
            for i in range(4):
                if self.q[i] <= x < self.q[i + 1]:
                    k = i
                    break
        # Increment positions above the cell.
        for i in range(k + 1, 5):
            self.npos[i] += 1
        for i in range(5):
            self.nprime[i] += self.dn[i]
        # Adjust interior markers.
        for i in range(1, 4):
            d = self.nprime[i] - self.npos[i]
            if (d >= 1 and self.npos[i + 1] - self.npos[i] > 1) or \
               (d <= -1 and self.npos[i - 1] - self.npos[i] < -1):
                sd = 1 if d >= 0 else -1
                qp = self._parabolic(i, sd)
                if self.q[i - 1] < qp < self.q[i + 1]:
                    self.q[i] = qp
                else:
                    self.q[i] = self._linear(i, sd)
                self.npos[i] += sd
        self.n += 1
        return self

    def _parabolic(self, i, d):
        n = self.npos
        q = self.q
        return q[i] + d / (n[i + 1] - n[i - 1]) * (
            (n[i] - n[i - 1] + d) * (q[i + 1] - q[i]) / (n[i + 1] - n[i])
            + (n[i + 1] - n[i] - d) * (q[i] - q[i - 1]) / (n[i] - n[i - 1]))

    def _linear(self, i, d):
        return self.q[i] + d * (self.q[i + d] - self.q[i]) / (self.npos[i + d] - self.npos[i])

    def value(self):
        """Current quantile estimate (exact order statistic while fewer than 5 seen)."""
        if self.n == 0:
            raise ValueError("no data")
        if self.n < 5:
            s = sorted(self.q)
            idx = min(int(self.p * self.n), self.n - 1)
            return s[idx]
        return self.q[2]


def _lcg_unit(seed):
    # Return a float in [0, 1) from the high bits (the low bits of this LCG are badly
    # non-random -- the lowest bit merely alternates -- so never use `state % m`).
    state = seed & 0x7FFFFFFF
    def rand():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return state / 0x80000000
    return rand


def reservoir_sample(stream, k, seed=1234567):
    """Vitter reservoir sample: ``k`` uniform items from a stream of unknown length.

    ``stream`` is any iterable. Returns a list of up to ``k`` items, each element of the
    stream equally likely to be included. One pass, O(k) memory. Deterministic for a
    fixed ``seed``.
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    rand = _lcg_unit(seed)
    reservoir = []
    for i, item in enumerate(stream):
        if i < k:
            reservoir.append(item)
        else:
            # Index from the high-bit float, not a low-bit modulo.
            j = int(rand() * (i + 1))
            if j < k:
                reservoir[j] = item
    return reservoir
