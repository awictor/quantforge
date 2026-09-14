"""Walker's alias method for O(1) categorical sampling.

Sampling from a fixed discrete distribution over ``n`` outcomes naively costs ``O(n)`` per
draw (a cumulative-probability search). Walker's alias method preprocesses the weights in
``O(n)`` into two tables so each subsequent draw is ``O(1)`` -- ideal when you draw many
times from the same distribution (Monte-Carlo scenario selection, weighted bootstrap).
Uses a deterministic LCG stream for reproducibility. Pure standard library.
"""


class AliasSampler:
    """O(1)-per-draw sampler for a fixed categorical distribution (Walker's alias method).

    Construct from a list of non-negative ``weights`` (need not sum to 1; normalized
    internally). ``sample()`` returns one index; ``sample_many(k)`` returns ``k`` indices.
    A ``seed`` fixes the deterministic random stream for reproducibility.
    """

    def __init__(self, weights, seed=1234567):
        n = len(weights)
        if n == 0:
            raise ValueError("need at least one weight")
        total = 0.0
        for w in weights:
            if w < 0:
                raise ValueError("weights must be non-negative")
            total += w
        if total <= 0.0:
            raise ValueError("weights must sum to a positive value")
        self.n = n
        self.prob = [0.0] * n
        self.alias = [0] * n
        # Scaled probabilities (mean 1).
        scaled = [w * n / total for w in weights]
        small = [i for i in range(n) if scaled[i] < 1.0]
        large = [i for i in range(n) if scaled[i] >= 1.0]
        while small and large:
            s = small.pop()
            l = large.pop()
            self.prob[s] = scaled[s]
            self.alias[s] = l
            scaled[l] = scaled[l] - (1.0 - scaled[s])
            if scaled[l] < 1.0:
                small.append(l)
            else:
                large.append(l)
        for i in large:
            self.prob[i] = 1.0
        for i in small:
            self.prob[i] = 1.0
        self._state = seed & 0x7FFFFFFF or 1

    def _rand(self):
        self._state = (1103515245 * self._state + 12345) & 0x7FFFFFFF
        return self._state / 0x80000000

    def sample(self):
        """Draw one index in ``[0, n)`` in ``O(1)``."""
        i = int(self._rand() * self.n)
        if i >= self.n:
            i = self.n - 1
        return i if self._rand() < self.prob[i] else self.alias[i]

    def sample_many(self, k):
        """Draw ``k`` indices as a list."""
        if k < 0:
            raise ValueError("k must be non-negative")
        return [self.sample() for _ in range(k)]
