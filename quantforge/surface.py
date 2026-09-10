"""Term-structure volatility surface: SVI fitted per expiry.

A single smile (see :mod:`quantforge.svi`) describes one expiry. A *surface*
stitches several expiries together and must be free of **calendar arbitrage**:
total implied variance ``w(k, t) = sigma^2 * t`` must be non-decreasing in
maturity at every fixed log-moneyness ``k``. If it ever falls as ``t`` grows,
a calendar spread locks in a riskless profit.

``VolSurface`` fits one raw-SVI slice per expiry, then interpolates total
variance linearly in ``t`` between slices (the standard no-arbitrage-friendly
scheme) and back out an implied vol at any ``(k, t)``. ``calendar_arbitrage``
reports any maturity pair whose variance curves cross.
"""

import math
from bisect import bisect_left
from dataclasses import dataclass
from typing import Sequence, List, Tuple

from .svi import SVIParams, calibrate_svi


@dataclass(frozen=True)
class SurfaceSlice:
    t: float                 # expiry in years
    params: SVIParams
    rmse: float


@dataclass(frozen=True)
class CalendarViolation:
    t_short: float
    t_long: float
    k: float
    w_short: float
    w_long: float           # < w_short => arbitrage


class VolSurface:
    """A term structure of SVI smiles with calendar-arbitrage diagnostics."""

    def __init__(self, slices: List[SurfaceSlice]):
        if not slices:
            raise ValueError("surface needs at least one slice")
        self.slices = sorted(slices, key=lambda s: s.t)
        self._ts = [s.t for s in self.slices]

    @classmethod
    def fit(cls, quotes, beta_ignored=None) -> "VolSurface":
        """Fit a surface from per-expiry quotes.

        Args:
            quotes: iterable of ``(t, ks, total_variances)`` triples, one per
                expiry, where ``ks`` are log-moneyness points and
                ``total_variances`` the observed w = sigma^2 * t.
        """
        slices = []
        for t, ks, tv in quotes:
            params, rmse = calibrate_svi(ks, tv)
            slices.append(SurfaceSlice(t=float(t), params=params, rmse=rmse))
        return cls(slices)

    def total_variance(self, k: float, t: float) -> float:
        """Interpolate total implied variance w(k, t) across the term structure.

        Linear in ``t`` between bracketing slices; flat-extrapolated in
        variance-per-year beyond the ends (constant forward variance).
        """
        ts = self._ts
        if t <= ts[0]:
            # Scale the first slice's variance down proportionally in t
            # (constant instantaneous variance back to zero).
            w0 = self.slices[0].params.total_variance(k)
            return w0 * (t / ts[0]) if ts[0] > 0 else w0
        if t >= ts[-1]:
            # Extrapolate at the last slice's forward variance rate.
            w_last = self.slices[-1].params.total_variance(k)
            if len(ts) >= 2:
                w_prev = self.slices[-2].params.total_variance(k)
                rate = (w_last - w_prev) / (ts[-1] - ts[-2])
            else:
                rate = w_last / ts[-1]
            return w_last + rate * (t - ts[-1])
        i = bisect_left(ts, t)
        t_lo, t_hi = ts[i - 1], ts[i]
        w_lo = self.slices[i - 1].params.total_variance(k)
        w_hi = self.slices[i].params.total_variance(k)
        frac = (t - t_lo) / (t_hi - t_lo)
        return w_lo + frac * (w_hi - w_lo)

    def implied_vol(self, k: float, t: float) -> float:
        """Black-Scholes implied vol at log-moneyness ``k`` and expiry ``t``."""
        if t <= 0:
            raise ValueError("t must be positive")
        w = self.total_variance(k, t)
        return math.sqrt(max(w, 0.0) / t)

    def calendar_arbitrage(self, ks: Sequence[float] = None,
                           tol: float = 1e-8) -> List[CalendarViolation]:
        """Check total variance is non-decreasing in ``t`` at each ``k``.

        Returns a list of violations (empty if the surface is calendar-arb
        free). Checks adjacent slice pairs on a grid of log-moneyness points.
        """
        if ks is None:
            ks = [-0.5, -0.25, -0.1, 0.0, 0.1, 0.25, 0.5]
        violations = []
        for a, b in zip(self.slices[:-1], self.slices[1:]):
            for k in ks:
                w_short = a.params.total_variance(k)
                w_long = b.params.total_variance(k)
                if w_long < w_short - tol:
                    violations.append(CalendarViolation(
                        t_short=a.t, t_long=b.t, k=k,
                        w_short=w_short, w_long=w_long))
        return violations

    def is_calendar_arbitrage_free(self, ks: Sequence[float] = None) -> bool:
        return not self.calendar_arbitrage(ks)

    def forward_variance(self, k: float, t1: float, t2: float) -> float:
        """Forward (instantaneous-average) variance between ``t1`` and ``t2``.

        The total variance is additive in time, so the variance realized over
        ``[t1, t2]`` at log-moneyness ``k`` is
        ``(w(k, t2) - w(k, t1)) / (t2 - t1)``. A negative result signals a
        calendar-arbitrage violation between those maturities.
        """
        if t2 <= t1:
            raise ValueError("require t2 > t1")
        w1 = self.total_variance(k, t1)
        w2 = self.total_variance(k, t2)
        return (w2 - w1) / (t2 - t1)

    def forward_vol(self, k: float, t1: float, t2: float) -> float:
        """Forward volatility between ``t1`` and ``t2`` = sqrt(forward variance).

        The vol of a forward-starting option that sets at ``t1`` and expires at
        ``t2``. Raises if the forward variance is negative (calendar arbitrage).
        """
        fv = self.forward_variance(k, t1, t2)
        if fv < 0:
            raise ValueError(
                f"negative forward variance ({fv:.4g}) between {t1} and {t2}; "
                "surface has calendar arbitrage there"
            )
        return math.sqrt(fv)
