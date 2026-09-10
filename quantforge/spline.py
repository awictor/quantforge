"""Natural cubic-spline interpolation (pure standard library).

A model-free alternative to a parametric smile (SVI/SABR): fit a natural cubic
spline through observed implied vols across strikes and read the vol at any
strike. The spline is C2-continuous with zero second derivative at the ends
(the "natural" boundary), which avoids the wild extrapolation a high-order
polynomial would give.

``CubicSpline`` is a general 1-D interpolator; ``SmileSpline`` wraps it for the
(strike -> implied vol) use case with flat extrapolation beyond the quoted
strikes (a common, safe smile convention).
"""

import bisect
from typing import Sequence, List


class CubicSpline:
    """Natural cubic spline through ``(xs, ys)`` with ``xs`` strictly increasing."""

    def __init__(self, xs: Sequence[float], ys: Sequence[float]):
        n = len(xs)
        if n < 3:
            raise ValueError("need at least three points for a cubic spline")
        if len(ys) != n:
            raise ValueError("xs and ys must be the same length")
        for i in range(n - 1):
            if xs[i + 1] <= xs[i]:
                raise ValueError("xs must be strictly increasing")
        self.xs = list(map(float, xs))
        self.ys = list(map(float, ys))
        self.m = self._second_derivatives()

    def _second_derivatives(self) -> List[float]:
        """Solve the tridiagonal system for the natural-spline moments m_i."""
        xs, ys = self.xs, self.ys
        n = len(xs)
        h = [xs[i + 1] - xs[i] for i in range(n - 1)]

        # Tridiagonal system A m = d for interior points (m_0 = m_{n-1} = 0).
        a = [0.0] * n  # sub-diagonal
        b = [1.0] * n  # diagonal
        c = [0.0] * n  # super-diagonal
        d = [0.0] * n  # rhs
        for i in range(1, n - 1):
            a[i] = h[i - 1]
            b[i] = 2.0 * (h[i - 1] + h[i])
            c[i] = h[i]
            d[i] = 6.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1])

        # Thomas algorithm.
        cp = [0.0] * n
        dp = [0.0] * n
        cp[0] = c[0] / b[0]
        dp[0] = d[0] / b[0]
        for i in range(1, n):
            denom = b[i] - a[i] * cp[i - 1]
            cp[i] = c[i] / denom if i < n - 1 else 0.0
            dp[i] = (d[i] - a[i] * dp[i - 1]) / denom
        m = [0.0] * n
        m[n - 1] = dp[n - 1]
        for i in range(n - 2, -1, -1):
            m[i] = dp[i] - cp[i] * m[i + 1]
        # Natural boundary conditions.
        m[0] = 0.0
        m[n - 1] = 0.0
        return m

    def __call__(self, x: float) -> float:
        """Evaluate the spline at ``x`` (clamped to the interpolation range)."""
        xs, ys, m = self.xs, self.ys, self.m
        if x <= xs[0]:
            return ys[0]
        if x >= xs[-1]:
            return ys[-1]
        i = bisect.bisect_right(xs, x) - 1
        i = max(0, min(i, len(xs) - 2))
        h = xs[i + 1] - xs[i]
        t = x - xs[i]
        # Standard cubic-spline segment formula.
        a_i = ys[i]
        b_i = (ys[i + 1] - ys[i]) / h - h * (2.0 * m[i] + m[i + 1]) / 6.0
        c_i = m[i] / 2.0
        d_i = (m[i + 1] - m[i]) / (6.0 * h)
        return a_i + b_i * t + c_i * t * t + d_i * t * t * t


class SmileSpline:
    """Strike -> implied-vol natural cubic spline with flat extrapolation."""

    def __init__(self, strikes: Sequence[float], vols: Sequence[float]):
        self._spline = CubicSpline(strikes, vols)
        self._lo_k, self._hi_k = strikes[0], strikes[-1]
        self._lo_v, self._hi_v = vols[0], vols[-1]

    def vol(self, K: float) -> float:
        """Implied vol at strike ``K``; flat beyond the quoted strike range."""
        if K <= self._lo_k:
            return self._lo_v
        if K >= self._hi_k:
            return self._hi_v
        return self._spline(K)
