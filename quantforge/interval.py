"""Interval arithmetic: rigorous enclosures under uncertainty.

An ``Interval`` represents a range ``[lo, hi]`` of possible values and propagates it
through arithmetic so the result is *guaranteed* to contain the true value -- the basis of
verified computation, tolerance analysis, and worst-case bounds. Operations follow the
standard interval rules (add endpoint-wise, multiply by taking min/max of the four
endpoint products, etc.). Pure standard library.

This uses ordinary floating point and does not perform outward rounding, so the enclosure
is mathematically correct but not a certified floating-point bound at the ULP level.
"""

import math


class Interval:
    """A closed interval ``[lo, hi]`` with guaranteed-enclosure arithmetic.

    Construct from two endpoints (order-normalized) or a single number (a degenerate
    point interval). Supports ``+ - * /`` with intervals and scalars, ``width``,
    ``midpoint``, ``contains``, ``intersect``, and monotone functions ``exp``/``log``/
    ``sqrt``/``__pow__`` (integer powers).
    """

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        if hi is None:
            hi = lo
        self.lo = float(min(lo, hi))
        self.hi = float(max(lo, hi))

    def __repr__(self):
        return "Interval(%g, %g)" % (self.lo, self.hi)

    def __eq__(self, other):
        return isinstance(other, Interval) and self.lo == other.lo and self.hi == other.hi

    @staticmethod
    def _coerce(v):
        return v if isinstance(v, Interval) else Interval(v)

    def __add__(self, other):
        o = self._coerce(other)
        return Interval(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __sub__(self, other):
        o = self._coerce(other)
        return Interval(self.lo - o.hi, self.hi - o.lo)

    def __rsub__(self, other):
        return self._coerce(other).__sub__(self)

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __mul__(self, other):
        o = self._coerce(other)
        products = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self._coerce(other)
        if o.lo <= 0.0 <= o.hi:
            raise ValueError("division by an interval containing zero")
        quotients = [self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi]
        return Interval(min(quotients), max(quotients))

    def __pow__(self, n):
        if not isinstance(n, int) or n < 0:
            raise ValueError("power must be a non-negative integer")
        if n == 0:
            return Interval(1.0, 1.0)
        lo_n = self.lo ** n
        hi_n = self.hi ** n
        if n % 2 == 1:
            return Interval(lo_n, hi_n)                 # odd power is monotone
        # even power: minimum is 0 if the interval straddles 0
        if self.lo <= 0.0 <= self.hi:
            return Interval(0.0, max(lo_n, hi_n))
        return Interval(min(lo_n, hi_n), max(lo_n, hi_n))

    def width(self):
        """Width ``hi - lo`` of the interval (the uncertainty)."""
        return self.hi - self.lo

    def midpoint(self):
        """Midpoint ``(lo + hi) / 2``."""
        return 0.5 * (self.lo + self.hi)

    def contains(self, value):
        """Whether a scalar ``value`` lies within ``[lo, hi]``."""
        return self.lo <= value <= self.hi

    def intersect(self, other):
        """Intersection with another interval, or ``None`` if they are disjoint."""
        lo = max(self.lo, other.lo)
        hi = min(self.hi, other.hi)
        return Interval(lo, hi) if lo <= hi else None

    def exp(self):
        """``exp`` of the interval (monotone increasing)."""
        return Interval(math.exp(self.lo), math.exp(self.hi))

    def log(self):
        """``log`` of a strictly positive interval (monotone increasing)."""
        if self.lo <= 0.0:
            raise ValueError("log requires a strictly positive interval")
        return Interval(math.log(self.lo), math.log(self.hi))

    def sqrt(self):
        """``sqrt`` of a non-negative interval (monotone increasing)."""
        if self.lo < 0.0:
            raise ValueError("sqrt requires a non-negative interval")
        return Interval(math.sqrt(self.lo), math.sqrt(self.hi))
