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

from .mathfns import norm_ppf, norm_pdf
from .bsm import delta as bs_delta, OptionType, price as bs_price
from .implied import implied_volatility
from .sabr import _solve3


def d1_d2(F, K, t, sigma):
    """Product ``d1 * d2`` of the Black d-terms on the forward ``F``."""
    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log(F / K) + 0.5 * vsqrt * vsqrt) / vsqrt
    return d1 * (d1 - vsqrt)


def _vega_vanna_volga(S, K, t, r, b, sigma):
    """Black-Scholes vega, vanna and volga at one strike (carry form).

    vanna = d vega / d spot, volga = d vega / d vol; all evaluated at ``sigma``.
    These three form the vanna-volga hedging basis.
    """
    vsqrt = sigma * math.sqrt(t)
    d1 = (math.log(S / K) + (b + 0.5 * sigma * sigma) * t) / vsqrt
    d2 = d1 - vsqrt
    carry = math.exp((b - r) * t)
    pdf = norm_pdf(d1)
    vega = S * carry * pdf * math.sqrt(t)
    vanna = -carry * pdf * d2 / sigma
    volga = vega * d1 * d2 / sigma
    return vega, vanna, volga


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

    def vol_at_delta(self, delta, is_call, tol=1e-10, max_iter=100):
        """Smile vol at a target (forward) Black-Scholes delta.

        The strike of a given delta depends on the vol at that strike, which the
        smile itself sets, so this solves the fixed point
        ``sigma = vol(strike_from_delta(sigma))`` by iteration seeded at the ATM
        vol. ``delta`` is the unsigned delta magnitude (e.g. ``0.25``); a put uses
        ``is_call=False``. At the pillar deltas it returns the pillar vols.
        """
        if not 0.0 < delta < 1.0:
            raise ValueError("delta magnitude must lie in (0, 1)")
        target = delta if is_call else -delta
        sigma = self.atm
        for _ in range(max_iter):
            K = _strike_from_delta(self.F, self.t, sigma, target, is_call)
            new = self.vol(K)
            if abs(new - sigma) < tol:
                return new
            sigma = new
        return sigma

    def risk_reversal(self, delta=0.25):
        """Smile-implied risk reversal at a given delta: ``sigma_call - sigma_put``.

        Recomputes the delta-consistent call and put vols from the fitted smile;
        at ``delta = 0.25`` this returns the input ``rr`` used to build the smile.
        """
        return (self.vol_at_delta(delta, True)
                - self.vol_at_delta(delta, False))

    def butterfly(self, delta=0.25):
        """Smile-implied butterfly at a given delta: ``(sigma_call + sigma_put)/2 - atm``.

        At ``delta = 0.25`` this returns the input ``bf`` used to build the smile.
        """
        return 0.5 * (self.vol_at_delta(delta, True)
                      + self.vol_at_delta(delta, False)) - self.atm

    def _cm_weights(self, K):
        """Vanna-volga replication weights x1, x2, x3 for strike ``K``.

        The weights make the three pillar options' vega, vanna and volga (at the
        flat ATM vol) replicate the target strike's -- Castagna-Mercurio's
        analytic result gives them in closed form from the pillar log-moneynesses.
        """
        k1, k2, k3 = self._ks
        x = math.log(K)
        y1, y2, y3 = math.log(k1), math.log(k2), math.log(k3)
        x1 = (x - y2) * (x - y3) / ((y1 - y2) * (y1 - y3))
        x2 = (x - y1) * (x - y3) / ((y2 - y1) * (y2 - y3))
        x3 = (x - y1) * (x - y2) / ((y3 - y1) * (y3 - y2))
        return x1, x2, x3

    def vol_cm(self, K, order=2):
        """Castagna-Mercurio vanna-volga implied-vol approximation at ``K``.

        ``order=1`` gives the first-order approximation -- a vega-weighted
        average of the pillar vols, ``sigma_atm + sum_i x_i (sigma_i -
        sigma_atm)`` with the vanna-volga replication weights. ``order=2`` adds
        the standard second-order correction

            sigma ~= sigma_atm + (-sigma_atm
                     + sqrt(sigma_atm^2 + d1 d2 (2 sigma_atm D1 + D2))) / (d1 d2),

        with ``D1`` the first-order excess and ``D2`` the pillar convexity term.
        Both are exact at the three pillars; the second order is the market
        standard.
        """
        k1, k2, k3 = self._ks
        s1, s2, s3 = self.sig[k1], self.sig[k2], self.sig[k3]
        s_atm = self.atm
        x1, x2, x3 = self._cm_weights(K)
        D1 = x1 * (s1 - s_atm) + x2 * (s2 - s_atm) + x3 * (s3 - s_atm)
        if order <= 1:
            return s_atm + D1

        vsqrt = s_atm * math.sqrt(self.t)
        d1 = (math.log(self.F / K) + 0.5 * vsqrt * vsqrt) / vsqrt
        d2 = d1 - vsqrt
        # Second-order convexity term: vanna-volga weighted squared vol excess.
        D2 = (x1 * d1_d2(self.F, k1, self.t, s_atm) * (s1 - s_atm) ** 2
              + x2 * d1_d2(self.F, k2, self.t, s_atm) * (s2 - s_atm) ** 2
              + x3 * d1_d2(self.F, k3, self.t, s_atm) * (s3 - s_atm) ** 2)
        disc = d1 * d2
        inner = s_atm * s_atm + disc * (2.0 * s_atm * D1 + D2)
        if inner < 0.0:
            return s_atm + D1        # fall back to first order if degenerate
        return s_atm + (-s_atm + math.sqrt(inner)) / disc if abs(disc) > 1e-12 \
            else s_atm + D1

    def price(self, K, r_dom=0.0, r_for=0.0, option_type=OptionType.CALL):
        """Vanna-volga option price via the Castagna-Mercurio correction.

        The exact second-order vanna-volga construction: start from the flat-ATM
        Black-Scholes price and add a linear combination of the three pillar
        options' *market-minus-ATM* price differences, weighted so the target
        option's vega, vanna and volga are matched by the hedging portfolio.
        This reprices the three market instruments exactly (unlike the quadratic
        vol interpolation in :meth:`vol`, which only matches the pillar vols) and
        is the standard FX smile pricer.

        ``r_dom``/``r_for`` are the domestic/foreign rates (carry ``b = r_dom -
        r_for``); pass the same pair used to build the forward.
        """
        b = r_dom - r_for
        s_atm = self.atm
        ks = self._ks

        # Target and pillar Greeks at the flat ATM vol (the hedging basis).
        vt, nt, ct = _vega_vanna_volga(self.S, K, self.t, r_dom, b, s_atm)
        rows = []  # columns: [vega, vanna, volga] of each pillar at ATM vol
        for ki in ks:
            rows.append(list(_vega_vanna_volga(self.S, ki, self.t, r_dom, b, s_atm)))

        # Solve for weights x so sum_i x_i * greeks(k_i) = greeks(K), i.e. the
        # pillar portfolio replicates the target's vega/vanna/volga.
        A = [[rows[i][j] for i in range(3)] for j in range(3)]  # 3x3, rows=greeks
        x = _solve3(A, [vt, nt, ct])
        if x is None:
            # Degenerate basis (e.g. coincident strikes): fall back to flat ATM.
            return bs_price(self.S, K, self.t, r_dom, s_atm, option_type, b=b)

        # Cost of the correction: each pillar's (market vol - ATM vol) price gap.
        correction = 0.0
        for i, ki in enumerate(ks):
            mkt = bs_price(self.S, ki, self.t, r_dom, self.sig[ki], option_type, b=b)
            atm = bs_price(self.S, ki, self.t, r_dom, s_atm, option_type, b=b)
            correction += x[i] * (mkt - atm)

        base = bs_price(self.S, K, self.t, r_dom, s_atm, option_type, b=b)
        return base + correction

    def vol_price_corrected(self, K, r_dom=0.0, r_for=0.0):
        """Implied vol of the Castagna-Mercurio vanna-volga price at ``K``.

        Inverts :meth:`price` back to a Black-Scholes vol, giving the smile that
        actually reprices the three market instruments. Slower than :meth:`vol`
        (one implied-vol solve per strike) but exact at the pillars in *price*,
        not merely in interpolated vol.
        """
        b = r_dom - r_for
        c = self.price(K, r_dom, r_for, OptionType.CALL)
        return implied_volatility(c, self.S, K, self.t, r_dom,
                                  OptionType.CALL, b=b)
