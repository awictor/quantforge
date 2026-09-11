"""Ho-Lee (1986) short-rate model: zero-coupon bonds.

Ho-Lee is the simplest no-arbitrage short-rate model: the rate has a
time-dependent drift and constant volatility, with no mean reversion:

    dr = theta(t) dt + sigma dW.

With a flat initial short rate ``r0`` (the constant-drift special case used
here, ``theta(t) = theta``), the zero-coupon bond has the affine closed form

    P(0, t) = exp(-r0 t - 0.5 theta t^2 + sigma^2 t^3 / 6),

so ``B(t) = t`` (the rate has unit sensitivity per year of maturity) and the
volatility contributes a ``+sigma^2 t^3 / 6`` convexity lift. Pure standard
library.
"""

import math

from .mathfns import norm_cdf


def holee_zero_coupon_bond(r0, t, theta, sigma):
    """Ho-Lee zero-coupon bond price P(0, t) with constant drift ``theta``.

    ``P = exp(-r0 t - 0.5 theta t^2 + sigma^2 t^3 / 6)``.
    """
    if t < 0:
        raise ValueError("t must be non-negative")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if t == 0:
        return 1.0
    return math.exp(-r0 * t - 0.5 * theta * t * t + sigma * sigma * t ** 3 / 6.0)


def holee_zero_coupon_yield(r0, t, theta, sigma):
    """Continuously-compounded yield of the Ho-Lee zero-coupon bond to ``t``.

    ``y(t) = r0 + 0.5 theta t - sigma^2 t^2 / 6`` (linear-in-t drift, quadratic
    convexity pull-down).
    """
    if t <= 0:
        raise ValueError("t must be positive")
    return r0 + 0.5 * theta * t - sigma * sigma * t * t / 6.0


def holee_expected_rate(r0, t, theta, sigma=0.0):
    """Expected Ho-Lee short rate ``E[r_t] = r0 + theta t``.

    With constant drift and no mean reversion the rate is
    ``r_t = r0 + theta t + sigma W_t``, so the mean drifts linearly at rate
    ``theta`` (``sigma`` does not enter the mean).
    """
    if t < 0:
        raise ValueError("t must be non-negative")
    return r0 + theta * t


def holee_rate_variance(t, sigma):
    """Variance of the Ho-Lee short rate ``Var[r_t] = sigma^2 t``.

    The rate is a drifted Brownian motion, so its variance grows linearly and
    without bound -- there is no stationary distribution (no mean reversion).
    """
    if t < 0:
        raise ValueError("t must be non-negative")
    return sigma * sigma * t


def holee_bond_greeks(r0, t, theta, sigma):
    """Exact rate sensitivities of a Ho-Lee zero-coupon bond.

    ``P = exp(-r0 t - ...)`` is linear in ``r0`` inside the exponent, so
    ``rho_r = dP/dr0 = -t P``, ``gamma_r = d2P/dr0^2 = t^2 P``, and the rate
    ``duration`` is exactly ``t`` (a Ho-Lee bond has duration equal to its
    maturity), with ``convexity = t^2``. Returns a dict with ``price``,
    ``rho_r``, ``gamma_r``, ``duration``, ``convexity``.
    """
    if t < 0:
        raise ValueError("t must be non-negative")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    price = holee_zero_coupon_bond(r0, t, theta, sigma)
    return {"price": price, "rho_r": -t * price, "gamma_r": t * t * price,
            "duration": t, "convexity": t * t}


def holee_bond_option(r0, t_option, t_bond, strike, theta, sigma, is_call=True):
    """European option on a Ho-Lee zero-coupon bond (exact Black-style).

    ``ln P(t_option, t_bond)`` is Gaussian, so the option is a Black formula on
    the forward bond ``P(0, t_bond) / P(0, t_option)`` with bond volatility
    ``sigma_p = sigma * (t_bond - t_option) * sqrt(t_option)`` (the Ho-Lee
    ``B(tau) = tau`` gives the linear maturity factor).
    """
    if not (0 < t_option < t_bond):
        raise ValueError("require 0 < t_option < t_bond")
    P_bond = holee_zero_coupon_bond(r0, t_bond, theta, sigma)
    P_opt = holee_zero_coupon_bond(r0, t_option, theta, sigma)
    sig_p = sigma * (t_bond - t_option) * math.sqrt(t_option)
    if sig_p < 1e-14:
        fwd = P_bond / P_opt
        payoff = max(fwd - strike, 0.0) if is_call else max(strike - fwd, 0.0)
        return P_opt * payoff
    d1 = (math.log(P_bond / (strike * P_opt)) + 0.5 * sig_p * sig_p) / sig_p
    d2 = d1 - sig_p
    if is_call:
        return P_bond * norm_cdf(d1) - strike * P_opt * norm_cdf(d2)
    return strike * P_opt * norm_cdf(-d2) - P_bond * norm_cdf(-d1)


def holee_coupon_bond_option(r0, t_option, cashflows, strike, theta, sigma,
                             is_call=True):
    """European option on a coupon bond under Ho-Lee (Jamshidian decomposition).

    ``cashflows`` is ``[(t_i, c_i), ...]`` with ``t_i > t_option``. The Ho-Lee
    bond is monotone decreasing in ``r0``, so Jamshidian's trick applies: solve
    for the critical rate ``r*`` where the coupon bond's value at expiry equals
    ``strike``, then sum the ``c_i``-weighted zero-coupon-bond options struck at
    ``K_i = P(t_option, t_i | r*)``. Exact.
    """
    cfs = sorted(cashflows)
    if not cfs or any(ti <= t_option for ti, _ in cfs):
        raise ValueError("all cashflow times must exceed t_option")

    def bond_value_at(r):
        return sum(c * holee_zero_coupon_bond(r, ti - t_option, theta, sigma)
                   for ti, c in cfs)

    lo, hi = -1.0, 1.0
    while bond_value_at(lo) < strike and lo > -50.0:
        lo -= 1.0
    while bond_value_at(hi) > strike and hi < 50.0:
        hi += 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if bond_value_at(mid) > strike:
            lo = mid
        else:
            hi = mid
    rstar = 0.5 * (lo + hi)

    total = 0.0
    for ti, c in cfs:
        Ki = holee_zero_coupon_bond(rstar, ti - t_option, theta, sigma)
        total += c * holee_bond_option(r0, t_option, ti, Ki, theta, sigma,
                                       is_call)
    return total


def holee_swaption(r0, expiry, pay_times, fixed_rate, theta, sigma, payer=True,
                   notional=1.0):
    """European swaption under Ho-Lee via the coupon-bond-option identity (exact).

    A payer swaption is a put on the fixed-leg coupon bond struck at the
    notional; a receiver is a call. Priced by :func:`holee_coupon_bond_option`.
    ``pay_times`` are the fixed-leg payment dates (all ``> expiry``); accruals
    are the gaps, the first measured from ``expiry``.
    """
    times = sorted(pay_times)
    if not times or any(t <= expiry for t in times):
        raise ValueError("all pay_times must exceed expiry")
    prev = expiry
    cfs = []
    for i, ti in enumerate(times):
        tau = ti - prev
        prev = ti
        c = fixed_rate * tau * notional
        if i == len(times) - 1:
            c += notional
        cfs.append((ti, c))
    is_call = not payer   # payer swaption = put on the coupon bond
    return holee_coupon_bond_option(r0, expiry, cfs, notional, theta, sigma,
                                    is_call)


def holee_caplet(r0, reset, pay, strike, theta, sigma, notional=1.0):
    """Caplet on ``[reset, pay]`` under Ho-Lee via the bond-put identity."""
    tau = pay - reset
    K_bond = 1.0 / (1.0 + strike * tau)
    put = holee_bond_option(r0, reset, pay, K_bond, theta, sigma, is_call=False)
    return notional * (1.0 + strike * tau) * put


def holee_floorlet(r0, reset, pay, strike, theta, sigma, notional=1.0):
    """Floorlet on ``[reset, pay]`` under Ho-Lee via the bond-call identity."""
    tau = pay - reset
    K_bond = 1.0 / (1.0 + strike * tau)
    call = holee_bond_option(r0, reset, pay, K_bond, theta, sigma, is_call=True)
    return notional * (1.0 + strike * tau) * call


def holee_cap(r0, dates, strike, theta, sigma, notional=1.0):
    """Ho-Lee cap: strip of caplets over successive ``dates`` (increasing times)."""
    ts = sorted(dates)
    if len(ts) < 2:
        raise ValueError("need >= 2 dates (one caplet)")
    return sum(holee_caplet(r0, ts[i - 1], ts[i], strike, theta, sigma, notional)
               for i in range(1, len(ts)))


def holee_floor(r0, dates, strike, theta, sigma, notional=1.0):
    """Ho-Lee floor: strip of floorlets over successive ``dates``."""
    ts = sorted(dates)
    if len(ts) < 2:
        raise ValueError("need >= 2 dates (one floorlet)")
    return sum(holee_floorlet(r0, ts[i - 1], ts[i], strike, theta, sigma,
                              notional) for i in range(1, len(ts)))
