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


def bootstrap_survival_curve(quote_maturities, quote_spreads, r, recovery=0.4,
                             freq=4, n_steps_per_year=100, tol=1e-10,
                             max_iter=100):
    """Bootstrap a piecewise-constant hazard curve from par CDS quotes.

    Given increasing ``quote_maturities`` and their par ``quote_spreads``, solve
    each tenor's forward hazard in turn (holding earlier segments fixed) so that
    the model par spread of :func:`cds_par_spread` reproduces the quote. Uses a
    bisection on the hazard, which is monotone in the par spread. Premium legs
    pay ``freq`` times a year; the protection-leg grid uses
    ``n_steps_per_year`` points per year. Returns the calibrated
    :class:`SurvivalCurve`.
    """
    n = len(quote_maturities)
    if len(quote_spreads) != n:
        raise ValueError("maturities and spreads must have equal length")
    if n == 0:
        raise ValueError("need at least one CDS quote")
    times = list(quote_maturities)
    hazards = []
    for i in range(n):
        T = times[i]
        pay = [k / freq for k in range(1, int(round(T * freq)) + 1)]
        n_steps = max(1, int(round(T * n_steps_per_year)))
        target = quote_spreads[i]

        def par_with(h):
            trial = SurvivalCurve(times[: i + 1], hazards + [h])
            return cds_par_spread(trial, pay, r, recovery, n_steps)

        # Par spread increases with the current-segment hazard: bracket + bisect.
        lo, hi = 1e-8, 1.0
        while par_with(hi) < target:
            hi *= 2.0
            if hi > 1e3:
                break
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            pm = par_with(mid)
            if abs(pm - target) < tol:
                lo = hi = mid
                break
            if pm < target:
                lo = mid
            else:
                hi = mid
        hazards.append(0.5 * (lo + hi))
    return SurvivalCurve(times, hazards)


def cds_greeks(curve: SurvivalCurve, spread, pay_times, r, recovery=0.4,
               n_steps=400, protection_buyer=True, bump=1e-4):
    """Risk sensitivities of a CDS mark-to-market by finite difference.

    Returns a dict with:

      * ``value``       -- the mark-to-market :func:`cds_value`;
      * ``credit01``    -- value change for a 1bp parallel bump of the hazard
        curve (credit spread risk);
      * ``ir01``        -- value change for a 1bp parallel bump of the discount
        rate ``r``;
      * ``recovery01``  -- value change for a 1-point (0.01) rise in recovery;
      * ``risky_annuity`` -- the survival-weighted premium annuity.

    A protection buyer gains when spreads widen (``credit01 > 0``) and loses as
    recovery rises. Bumps are one-sided by ``bump`` (hazard/rate) or 0.01
    (recovery).
    """
    base = cds_value(curve, spread, pay_times, r, recovery, n_steps,
                     protection_buyer)
    bumped_curve = SurvivalCurve(curve.times, [h + bump for h in curve.hazards])
    credit01 = (cds_value(bumped_curve, spread, pay_times, r, recovery, n_steps,
                          protection_buyer) - base)
    ir01 = (cds_value(curve, spread, pay_times, r + bump, recovery, n_steps,
                      protection_buyer) - base)
    rec01 = (cds_value(curve, spread, pay_times, r, recovery + 0.01, n_steps,
                       protection_buyer) - base)
    ann = risky_annuity(curve, pay_times, r)
    return {"value": base, "credit01": credit01, "ir01": ir01,
            "recovery01": rec01, "risky_annuity": ann}


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
