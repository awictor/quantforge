"""Vanna-Volga FX smile construction from three market quotes.

FX options markets quote the smile with three numbers per expiry: the
at-the-money vol (``atm``), the 25-delta risk reversal (``rr``, the call-minus-
put vol skew), and the 25-delta butterfly (``bf``, the smile convexity). From
these the three pillar vols are recovered:

    sigma_ATM = atm
    sigma_25C = atm + bf + rr/2      (25-delta call)
    sigma_25P = atm + bf - rr/2      (25-delta put)

The vanna-volga method then interpolates/extrapolates the vol at any strike as
a second-order correction to the ATM vol, using the three pillars. This module
builds the pillar strikes from their deltas and returns a callable smile.

Pure standard library; uses the Black-Scholes delta and the accurate inverse
normal CDF from :mod:`quantforge.mathfns`.
"""

import math

from .mathfns import norm_ppf
from .bsm import delta as bs_delta, OptionType


def pillar_vols(atm, rr, bf):
    """Return (sigma_25put, sigma_atm, sigma_25call) from ATM / RR / BF quotes."""
    return (atm + bf - 0.5 * rr, atm, atm + bf + 0.5 * rr)


def _strike_from_delta(F, t, sigma, target_delta, is_call):
    """Strike with the given (forward) delta under Black-Scholes.

    Uses the forward-delta convention: delta_call = N(d1), so
    K = F * exp(-sigma*sqrt(t)*N^{-1}(delta) + 0.5 sigma^2 t) for a call.
    """
    vsqrt = sigma * math.sqrt(t)
    if is_call:
        d1 = norm_ppf(target_delta)
    else:
        d1 = norm_ppf(1.0 + target_delta)   # put delta is negative
    return F * math.exp(-d1 * vsqrt + 0.5 * sigma * sigma * t)


class VannaVolgaSmile:
    """A vanna-volga FX smile built from ATM / RR / BF at one expiry."""

    def __init__(self, S, t, r_dom, r_for, atm, rr, bf, call_delta=0.25):
        self.S = S
        self.t = t
        self.F = S * math.exp((r_dom - r_for) * t)
        self.atm = atm
        s_p, s_a, s_c = pillar_vols(atm, rr, bf)
        # Pillar strikes from their 25-delta definitions.
        self.K_atm = self.F * math.exp(0.5 * s_a * s_a * t)   # delta-neutral ATM
        self.K_c = _strike_from_delta(self.F, t, s_c, call_delta, True)
        self.K_p = _strike_from_delta(self.F, t, s_p, -call_delta, False)
        self.sig = {self.K_p: s_p, self.K_atm: s_a, self.K_c: s_c}
        self._ks = sorted(self.sig)

    def vol(self, K):
        """Vanna-volga smile vol at strike ``K`` (second-order interpolation).

        Uses the standard second-order vanna-volga approximation: the vol is the
        ATM vol plus a quadratic-in-log-strike correction fit through the three
        pillars. Exact at the three pillar strikes.
        """
        k1, k2, k3 = self._ks
        s1, s2, s3 = self.sig[k1], self.sig[k2], self.sig[k3]
        x = math.log(K)
        x1, x2, x3 = math.log(k1), math.log(k2), math.log(k3)
        # Lagrange quadratic through (x_i, s_i) -- the smile in log-strike.
        L1 = (x - x2) * (x - x3) / ((x1 - x2) * (x1 - x3))
        L2 = (x - x1) * (x - x3) / ((x2 - x1) * (x2 - x3))
        L3 = (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2))
        return s1 * L1 + s2 * L2 + s3 * L3

    def pillars(self):
        """Return the three (strike, vol) pillar points, sorted by strike."""
        return [(k, self.sig[k]) for k in self._ks]
