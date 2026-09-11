"""Reduced-form credit: hazard-rate survival curve and CDS pricing.

A default time ``tau`` is modelled by a piecewise-constant hazard rate ``h``
(intensity). The survival probability to ``t`` is ``Q(t) = exp(-integral_0^t h)``
and the default density is ``h(t) Q(t)``. A credit default swap exchanges a
periodic premium (paid while the name survives) for a contingent
``(1 - recovery)`` payment at default:

    premium leg PV = spread * sum_i tau_i DF(t_i) Q(t_i)          (+ accrual)
    protection leg PV = (1 - R) * integral DF(t) (-dQ)            (approx by grid)
    par spread = protection PV / risky annuity.

Discounting uses a flat continuously-compounded rate ``r`` (or a supplied
discount function). Pure standard library.
"""

import math
from bisect import bisect_right
from typing import Sequence


class SurvivalCurve:
    """Piecewise-constant hazard-rate survival curve.

    Built from pillar times and the *forward* hazard rate on each segment
    ``[t_{i-1}, t_i]``. ``survival(t)`` returns ``Q(t) = exp(-integral h)`` and
    ``default_density(t)`` returns ``h(t) Q(t)``.
    """

    def __init__(self, times: Sequence[float], hazards: Sequence[float]):
        if len(times) != len(hazards):
            raise ValueError("times and hazards must have equal length")
        if not times:
            raise ValueError("need at least one pillar")
        prev = 0.0
        for t in times:
            if t <= prev:
                raise ValueError("times must be strictly increasing and positive")
            prev = t
        if any(h < 0 for h in hazards):
            raise ValueError("hazard rates must be non-negative")
        self.times = list(times)
        self.hazards = list(hazards)
        # Cumulative hazard at each pillar for O(log n) survival lookups.
        self._cum = []
        acc = 0.0
        prev = 0.0
        for t, h in zip(self.times, self.hazards):
            acc += h * (t - prev)
            self._cum.append(acc)
            prev = t

    def _cumulative_hazard(self, t):
        if t <= 0.0:
            return 0.0
        i = bisect_right(self.times, t)
        if i == 0:
            return self.hazards[0] * t
        base = self._cum[i - 1]
        if i < len(self.times):
            return base + self.hazards[i] * (t - self.times[i - 1])
        # Beyond the last pillar: extend with the final hazard.
        return base + self.hazards[-1] * (t - self.times[-1])

    def survival(self, t):
        """Survival probability ``Q(t) = P(tau > t)``."""
        return math.exp(-self._cumulative_hazard(t))

    def hazard(self, t):
        """Forward hazard rate applying at time ``t``."""
        if t <= self.times[0]:
            return self.hazards[0]
        i = bisect_right(self.times, t)
        return self.hazards[min(i, len(self.hazards) - 1)]

    def default_density(self, t):
        """Unconditional default density ``h(t) Q(t)``."""
        return self.hazard(t) * self.survival(t)


def _disc_fn(r):
    return r if callable(r) else (lambda t: math.exp(-r * t))


def risky_annuity(curve: SurvivalCurve, pay_times, r, accrual=None):
    """Risky (survival-weighted) annuity ``sum_i tau_i DF(t_i) Q(t_i)``.

    ``pay_times`` are the premium payment dates; ``accrual`` is the per-period
    year fractions (defaults to the gaps between pay times, starting from 0).
    Discount by flat rate ``r`` or a supplied ``r(t)`` function.
    """
    df = _disc_fn(r)
    a = 0.0
    prev = 0.0
    for k, t in enumerate(pay_times):
        tau = (t - prev) if accrual is None else accrual[k]
        a += tau * df(t) * curve.survival(t)
        prev = t
    return a


def cds_protection_leg(curve: SurvivalCurve, maturity, r, recovery=0.4,
                       n_steps=400):
    """PV of the CDS protection (default) leg, ``(1-R) integral DF(t) (-dQ)``.

    Numerically integrates the loss payment over ``[0, maturity]`` on a uniform
    grid, paying ``(1 - recovery)`` at the (grid-approximated) default time.
    """
    df = _disc_fn(r)
    dt = maturity / n_steps
    pv = 0.0
    q_prev = curve.survival(0.0)
    for i in range(1, n_steps + 1):
        t = i * dt
        q = curve.survival(t)
        # Default probability in (t-dt, t]; pay at the interval midpoint.
        pv += df(t - 0.5 * dt) * (q_prev - q)
        q_prev = q
    return (1.0 - recovery) * pv


def cds_premium_leg(curve: SurvivalCurve, spread, pay_times, r, accrual=None):
    """PV of the CDS premium leg at a given ``spread`` (annualized)."""
    return spread * risky_annuity(curve, pay_times, r, accrual)


def cds_par_spread(curve: SurvivalCurve, pay_times, r, recovery=0.4,
                   n_steps=400):
    """Fair (par) CDS spread: protection-leg PV divided by the risky annuity."""
    maturity = pay_times[-1]
    prot = cds_protection_leg(curve, maturity, r, recovery, n_steps)
    ann = risky_annuity(curve, pay_times, r)
    if ann <= 0:
        raise ValueError("risky annuity must be positive")
    return prot / ann


def cds_value(curve: SurvivalCurve, spread, pay_times, r, recovery=0.4,
              n_steps=400, protection_buyer=True):
    """Mark-to-market value of a CDS at a contractual ``spread``.

    Protection buyer is long the protection leg and short the premium leg:
    ``V = protection - spread * annuity``. Positive when the par spread has
    widened beyond the contractual spread.
    """
    maturity = pay_times[-1]
    prot = cds_protection_leg(curve, maturity, r, recovery, n_steps)
    prem = cds_premium_leg(curve, spread, pay_times, r)
    v = prot - prem
    return v if protection_buyer else -v
