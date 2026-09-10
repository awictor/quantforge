"""Discount curve: interpolation and a par-swap bootstrap.

A discount curve maps maturity ``T`` to the price today of $1 paid at ``T``,
``P(0, T) = DF(T)``. This module builds one from pillar discount factors (or
continuously-compounded zero rates) with log-linear interpolation in ``T`` --
equivalent to piecewise-constant instantaneous forward rates, the market-
standard no-arbitrage-friendly scheme -- and bootstraps pillar factors from a
set of par (fair) swap rates. The resulting callable ``DF(T)`` is exactly what
the Gaussian short-rate models (``cheyette``, ``g2pp``) consume as their initial
curve, so those models reprice the curve's own bonds exactly by construction.
Pure standard library.
"""

import math
from bisect import bisect_left


class DiscountCurve:
    """Log-linear discount curve from pillar ``(T, DF)`` points."""

    def __init__(self, times, dfs):
        pts = sorted(zip((float(t) for t in times), (float(d) for d in dfs)))
        self.times = [t for t, _ in pts]
        self.dfs = [d for _, d in pts]
        if not self.times or self.times[0] < 0:
            raise ValueError("need non-negative maturities")
        if any(d <= 0 for d in self.dfs):
            raise ValueError("discount factors must be positive")
        # Prepend T=0, DF=1 if absent so short-end interpolation is anchored.
        if self.times[0] > 1e-12:
            self.times.insert(0, 0.0)
            self.dfs.insert(0, 1.0)
        self._logdf = [math.log(d) for d in self.dfs]

    @classmethod
    def from_zero_rates(cls, times, zero_rates):
        """Build from continuously-compounded zero rates: ``DF = e^{-z T}``."""
        dfs = [math.exp(-z * t) for t, z in zip(times, zero_rates)]
        return cls(times, dfs)

    def df(self, T):
        """Discount factor ``P(0, T)`` by log-linear interpolation in ``T``."""
        T = float(T)
        if T <= 0.0:
            return 1.0
        ts = self.times
        if T >= ts[-1]:
            # Flat-forward extrapolation beyond the last pillar.
            f = (self._logdf[-1] - self._logdf[-2]) / (ts[-1] - ts[-2])
            return math.exp(self._logdf[-1] + f * (T - ts[-1]))
        i = bisect_left(ts, T)
        if ts[i] == T:
            return self.dfs[i]
        t0, t1 = ts[i - 1], ts[i]
        l0, l1 = self._logdf[i - 1], self._logdf[i]
        return math.exp(l0 + (l1 - l0) * (T - t0) / (t1 - t0))

    def __call__(self, T):
        return self.df(T)

    def zero_rate(self, T):
        """Continuously-compounded zero rate ``-ln DF(T) / T``."""
        if T <= 0.0:
            return 0.0
        return -math.log(self.df(T)) / T

    def forward_rate(self, T1, T2):
        """Continuously-compounded forward rate over ``[T1, T2]``."""
        if T2 <= T1:
            raise ValueError("need T2 > T1")
        return (math.log(self.df(T1)) - math.log(self.df(T2))) / (T2 - T1)

    def par_swap_rate(self, pay_times):
        """Par (fair fixed) rate of a swap with the given annual pay schedule.

        ``par = (1 - DF(T_n)) / sum_i tau_i DF(T_i)``, with the float leg valued
        as ``1 - DF(T_n)`` (unit notional, spot start).
        """
        annuity = 0.0
        prev = 0.0
        for Ti in pay_times:
            annuity += (Ti - prev) * self.df(Ti)
            prev = Ti
        return (1.0 - self.df(pay_times[-1])) / annuity


def bootstrap_from_swaps(swap_maturities, par_rates, freq=1.0):
    """Bootstrap a :class:`DiscountCurve` from par swap rates.

    Args:
        swap_maturities: increasing swap tenors in years (each an integer number
            of ``1/freq``-year periods).
        par_rates: the fair fixed rate for each tenor.
        freq: fixed-leg payments per year (1 = annual).

    Solves pillar by pillar: with all shorter discount factors known, each new
    par-rate equation is linear in the final ``DF(T_n)``. Returns the curve whose
    par-swap rates reproduce the inputs.
    """
    if len(swap_maturities) != len(par_rates):
        raise ValueError("maturities and par_rates must match in length")
    tau = 1.0 / freq
    times = []
    dfs = []

    def df_at(T):
        # Interpolate on the pillars fixed so far (log-linear); anchor DF(0)=1.
        if T <= 1e-12:
            return 1.0
        if not times:
            return 1.0
        if T <= times[0]:
            l = math.log(dfs[0]) * (T / times[0])
            return math.exp(l)
        if T >= times[-1]:
            return dfs[-1]
        i = bisect_left(times, T)
        if times[i] == T:
            return dfs[i]
        t0, t1 = times[i - 1], times[i]
        l0, l1 = math.log(dfs[i - 1]), math.log(dfs[i])
        return math.exp(l0 + (l1 - l0) * (T - t0) / (t1 - t0))

    for mat, par in zip(swap_maturities, par_rates):
        n = int(round(mat * freq))
        pay = [(k + 1) * tau for k in range(n)]
        # annuity over all-but-last period (known), plus the unknown last DF.
        known_ann = 0.0
        for Ti in pay[:-1]:
            known_ann += tau * df_at(Ti)
        Tn = pay[-1]
        # par = (1 - DF_n) / (known_ann + tau DF_n)  ->  solve for DF_n.
        DF_n = (1.0 - par * known_ann) / (1.0 + par * tau)
        times.append(Tn)
        dfs.append(DF_n)

    return DiscountCurve(times, dfs)
