"""Range-accrual note (closed form).

A range-accrual note pays a coupon proportional to the fraction of observation
dates on which a reference index stays inside a range ``[L, U]``. Over ``m``
equally weighted observation dates ``t_1 < ... < t_m <= T`` the accrued coupon is

    coupon * (1/m) * sum_i 1{ L <= S_{t_i} <= U }

paid at maturity ``T``. Under geometric Brownian motion each indicator has
risk-neutral probability

    P(L <= S_t <= U) = N(d(L)) - N(d(U)),   d(x) = [ln(S/x) + (b - sigma^2/2) t] / (sigma sqrt(t))

(the probability the *terminal* price sits in the band, using the real drift ``b``
rather than the risk-neutral pricing drift -- an accrual counts physical time in
the band). The present value is the discounted expected coupon, a closed-form sum
of these range probabilities. Pure standard library.
"""

import math

from .mathfns import norm_cdf


def _in_range_prob(S, L, U, t, b, sigma):
    """Risk-neutral probability that ``S_t`` lands in ``[L, U]`` under GBM."""
    vt = sigma * math.sqrt(t)
    drift = (b - 0.5 * sigma * sigma) * t
    # P(S_t <= x) = N( (ln(x/S) - drift) / vt ); band prob is the difference.
    dU = (math.log(U / S) - drift) / vt
    dL = (math.log(L / S) - drift) / vt
    return norm_cdf(dU) - norm_cdf(dL)


def range_accrual_note(S, L, U, t, r, sigma, coupon, observations, b=None, q=0.0,
                       notional=1.0):
    """Present value of a range-accrual note's coupon leg.

    Parameters
    ----------
    S, L, U : float
        Spot and the lower/upper edges of the accrual band, ``0 < L < U``.
    t : float
        Maturity in years; the coupon is paid at ``t``.
    r, sigma : float
        Risk-free rate and volatility.
    coupon : float
        Full coupon rate earned if the index is in range on every observation.
    observations : int
        Number of equally spaced observation dates in ``(0, t]``. Date ``i`` of
        ``m`` falls at ``t_i = t * i / m``.
    b : float, optional
        Cost of carry / index drift. Defaults to ``r - q``.
    q : float
        Dividend yield, used only when ``b`` is not given.
    notional : float
        Note notional.

    Returns
    -------
    float
        Discounted expected coupon. Non-negative, rises with a wider band, and
        approaches ``notional * coupon * e^{-r t}`` as the band widens to cover
        the whole positive axis.
    """
    if S <= 0 or t <= 0 or sigma <= 0:
        raise ValueError("S, t, sigma must be positive")
    if not (0.0 < L < U):
        raise ValueError("require 0 < L < U")
    if observations < 1:
        raise ValueError("observations must be >= 1")
    if b is None:
        b = r

    m = int(observations)
    frac = 0.0
    for i in range(1, m + 1):
        ti = t * i / m
        frac += _in_range_prob(S, L, U, ti, b, sigma)
    frac /= m
    return notional * coupon * frac * math.exp(-r * t)
