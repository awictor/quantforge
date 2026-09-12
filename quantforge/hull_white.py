"""Hull-White (extended Vasicek) one-factor short-rate model.

Hull-White adds a time-dependent drift to Vasicek,

    dr = (theta(t) - a r) dt + sigma dW,

so the model reprices *any* initial zero curve exactly while keeping analytic bond
prices. The fitted bond price is

    P(t, T) = A(t, T) exp(-B(t, T) r_t),   B = (1 - e^{-a(T-t)}) / a,

with ``A`` chosen from the market discount factors and the forward curve so that
``P(0, T)`` equals the input. This module works from an initial discount-factor
curve (a callable ``P0(T)``) and returns fitted zero prices and the ``B`` factor.
Pure standard library.
"""

import math


def hw_B(a, tau):
    """Hull-White ``B(t, T) = (1 - e^{-a*tau}) / a`` for ``tau = T - t``."""
    if a == 0.0:
        return tau
    return (1.0 - math.exp(-a * tau)) / a


def hw_zero_from_curve(P0, r0, a, sigma, t, T, f0=None, eps=1e-5):
    """Hull-White bond price ``P(t, T)`` fitted to an initial curve ``P0``.

    Parameters
    ----------
    P0 : callable
        Initial discount factor ``P0(T)`` observed today (``P0(0) = 1``).
    r0 : float
        Current short rate; for consistency it should equal the initial instant
        forward ``f(0,0)``.
    a, sigma : float
        Mean reversion and volatility (``a`` may be 0 for the Ho-Lee limit).
    t, T : float
        Valuation and maturity times, ``0 <= t <= T``.
    f0 : callable, optional
        Initial instantaneous forward ``f(0, t)``; defaults to a finite-difference
        of ``-ln P0``.
    eps : float
        Step for the forward finite difference.

    Returns
    -------
    float
        The fitted ``P(t, T)``. At ``t = 0`` it reproduces ``P0(T)`` exactly.
    """
    if sigma < 0 or t < 0 or T < t:
        raise ValueError("require sigma >= 0 and 0 <= t <= T")
    if f0 is None:
        def f0(u):
            pu = P0(max(u, 0.0))
            pv = P0(u + eps)
            return -(math.log(pv) - math.log(pu)) / eps

    B = hw_B(a, T - t)
    # Analytic A(t,T) that matches the initial curve (Hull-White 1994).
    lnA = (math.log(P0(T) / P0(t)) + B * f0(t)
           - (sigma * sigma) / (4.0 * a) * B * B * (1.0 - math.exp(-2.0 * a * t))
           if a != 0.0 else
           math.log(P0(T) / P0(t)) + B * f0(t)
           - 0.5 * sigma * sigma * t * B * B)
    return math.exp(lnA - B * r0)


def hw_bond_option(P0, a, sigma, t_option, t_bond, strike, is_call=True):
    """European option on a zero-coupon bond under Hull-White (analytic).

    Prices an option expiring at ``t_option`` on a zero maturing at ``t_bond``,
    struck at ``strike``, using the initial discount curve ``P0``. The forward
    bond price is log-normal with volatility

        sigma_P = sigma * B(a, t_bond - t_option) * sqrt((1 - e^{-2 a t_option}) / (2 a)),

    giving a Black-style formula in the discount factors ``P0(t_bond)`` and
    ``P0(t_option)`` (Jamshidian / Hull-White). Returns the option value today.
    Requires ``0 < t_option < t_bond``.
    """
    import math as _m
    from .mathfns import norm_cdf as _N

    if sigma < 0 or not (0.0 < t_option < t_bond):
        raise ValueError("require sigma >= 0 and 0 < t_option < t_bond")
    pT = P0(t_bond)
    pS = P0(t_option)
    if a != 0.0:
        var = (sigma * sigma) * (1.0 - _m.exp(-2.0 * a * t_option)) / (2.0 * a)
    else:
        var = sigma * sigma * t_option
    sig_p = hw_B(a, t_bond - t_option) * _m.sqrt(var)
    if sig_p <= 0.0:
        # Degenerate: intrinsic on the forward.
        fwd = pT / pS
        payoff = max(fwd - strike, 0.0) if is_call else max(strike - fwd, 0.0)
        return pS * payoff
    h = _m.log(pT / (pS * strike)) / sig_p + 0.5 * sig_p
    if is_call:
        return pT * _N(h) - strike * pS * _N(h - sig_p)
    return strike * pS * _N(-h + sig_p) - pT * _N(-h)


def hw_caplet(P0, a, sigma, reset, pay, strike, notional=1.0):
    """Hull-White caplet: an option on the simple forward rate over ``[reset, pay]``.

    A caplet paying ``notional * tau * max(L - strike, 0)`` at ``pay`` (where ``L``
    is the simple rate set at ``reset`` for accrual ``tau = pay - reset``) equals
    ``notional * (1 + strike*tau)`` puts on the ``pay``-zero struck at
    ``1/(1 + strike*tau)``, expiring at ``reset`` -- the standard bond-option
    representation. Priced analytically off the initial curve.
    """
    tau = pay - reset
    if tau <= 0:
        raise ValueError("require pay > reset")
    k = 1.0 / (1.0 + strike * tau)
    put = hw_bond_option(P0, a, sigma, reset, pay, k, is_call=False)
    return notional * (1.0 + strike * tau) * put


def hw_floorlet(P0, a, sigma, reset, pay, strike, notional=1.0):
    """Hull-White floorlet: ``notional * (1 + strike*tau)`` calls on the pay-zero."""
    tau = pay - reset
    if tau <= 0:
        raise ValueError("require pay > reset")
    k = 1.0 / (1.0 + strike * tau)
    call = hw_bond_option(P0, a, sigma, reset, pay, k, is_call=True)
    return notional * (1.0 + strike * tau) * call


def hw_cap(P0, a, sigma, dates, strike, notional=1.0):
    """Hull-White cap: sum of caplets over consecutive ``dates`` (reset, pay pairs).

    ``dates`` is the schedule ``[t_0, t_1, ..., t_n]``; caplet ``i`` covers
    ``[t_i, t_{i+1}]``. A floor is the analogous sum of floorlets. By put-call
    parity ``cap - floor`` equals the value of the fixed-vs-float swap.
    """
    if len(dates) < 2:
        raise ValueError("need at least two schedule dates")
    return sum(hw_caplet(P0, a, sigma, dates[i], dates[i + 1], strike, notional)
               for i in range(len(dates) - 1))


def hw_floor(P0, a, sigma, dates, strike, notional=1.0):
    """Hull-White floor: sum of floorlets over consecutive ``dates``."""
    if len(dates) < 2:
        raise ValueError("need at least two schedule dates")
    return sum(hw_floorlet(P0, a, sigma, dates[i], dates[i + 1], strike, notional)
               for i in range(len(dates) - 1))


def hw_swaption(P0, r0, a, sigma, expiry, pay_times, fixed_rate,
                is_payer=True, notional=1.0):
    """European swaption under Hull-White via Jamshidian decomposition.

    At ``expiry`` the holder enters a swap paying (payer) or receiving (receiver)
    ``fixed_rate`` on ``pay_times``. Jamshidian's trick: the underlying coupon bond
    is monotone in the short rate, so find the critical rate ``r*`` at which its
    value equals par, then the swaption is a portfolio of options on each
    zero-coupon cashflow struck at that cashflow's value at ``r*``. Priced off the
    initial curve; ``r0`` is only used to seed the ``r*`` search.

    A payer swaption is a put on the coupon bond (a portfolio of zero puts); a
    receiver is the corresponding call portfolio.
    """
    import math as _m
    if not pay_times or any(pay_times[i] <= expiry for i in range(len(pay_times))):
        raise ValueError("all pay_times must be after expiry")
    if sigma < 0:
        raise ValueError("sigma must be non-negative")

    # Cashflows of the fixed leg (coupon = fixed_rate * accrual), plus notional.
    prev = expiry
    cfs = []
    for pt in pay_times:
        tau = pt - prev
        cfs.append((pt, fixed_rate * tau))
        prev = pt
    cfs[-1] = (cfs[-1][0], cfs[-1][1] + 1.0)     # principal at the end

    # Bond value at expiry as a function of the short rate r (HW: P(expiry,T)
    # given r), using the analytic reconstitution around the fitted curve.
    def bond_at(r):
        total = 0.0
        for (T, c) in cfs:
            # P(expiry, T | r) = A(expiry,T) exp(-B r), matched to the curve.
            p = hw_zero_from_curve(P0, r, a, sigma, expiry, T)
            total += c * p
        return total

    # Critical r* where the coupon bond equals par (1.0). Bond is decreasing in r.
    lo, hi = -0.5, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if bond_at(mid) > 1.0:
            lo = mid
        else:
            hi = mid
    r_star = 0.5 * (lo + hi)

    # Strike for each zero option = its value at r*; sum the bond options.
    total = 0.0
    for (T, c) in cfs:
        k = hw_zero_from_curve(P0, r_star, a, sigma, expiry, T)
        # Payer = put on the bond = puts on each zero; receiver = calls.
        opt = hw_bond_option(P0, a, sigma, expiry, T, k, is_call=not is_payer)
        total += c * opt
    return notional * total
