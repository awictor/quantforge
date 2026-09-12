"""Structured notes: principal-protected notes and reverse convertibles.

Retail structured products decompose into a bond plus an option position:

    Principal-protected note = zero-coupon bond (returns principal at maturity)
                               + participation * call on the underlying
    Reverse convertible      = zero-coupon bond + coupon - short put

The PPN guarantees the principal (its value never falls below the discounted
principal) while giving equity upside through the call; the reverse convertible
pays an enhanced coupon funded by selling downside via the put. This module
values both from their components on top of :mod:`quantforge.bsm`. Pure standard
library.
"""

import math

from .bsm import call_price, put_price


def note_zero_coupon_bond(principal, r, maturity):
    """Present value of a zero-coupon bond: ``principal * e^{-r T}``."""
    if principal < 0:
        raise ValueError("principal must be non-negative")
    if maturity < 0:
        raise ValueError("maturity must be non-negative")
    return principal * math.exp(-r * maturity)


def principal_protected_note(S, K, maturity, r, sigma, principal,
                             participation=1.0, n_shares=None, b=None):
    """Value of a principal-protected note = ZC bond + participation * call.

    Guarantees the ``principal`` at maturity (the ZC-bond leg) and adds
    ``participation`` times a call on the underlying for upside. ``n_shares``
    scales the call to the note's notional (defaults to ``principal / S``, i.e. the
    note buys as many shares as the principal affords at inception). The value is
    always at least the discounted principal (the call leg is non-negative).
    """
    if participation < 0:
        raise ValueError("participation must be non-negative")
    if n_shares is None:
        if S <= 0:
            raise ValueError("S must be positive to default n_shares")
        n_shares = principal / S
    bond = note_zero_coupon_bond(principal, r, maturity)
    call = n_shares * call_price(S, K, maturity, r, sigma, b)
    return bond + participation * call


def reverse_convertible(S, K, maturity, r, sigma, principal, coupon_rate,
                        n_shares=None, b=None):
    """Value of a reverse convertible = ZC bond + coupon PV - short put.

    The investor receives an enhanced coupon and the principal, but is short a put
    struck at ``K`` (delivering shares if the underlying falls). Value is
    ``PV(principal + coupon) - n_shares * put``. ``n_shares`` defaults to
    ``principal / K`` (the put covers the principal at the strike). Below the
    plain bond-plus-coupon value because of the short put.
    """
    if coupon_rate < 0:
        raise ValueError("coupon_rate must be non-negative")
    if n_shares is None:
        if K <= 0:
            raise ValueError("K must be positive to default n_shares")
        n_shares = principal / K
    disc = math.exp(-r * maturity)
    bond_plus_coupon = principal * (1.0 + coupon_rate * maturity) * disc
    put = n_shares * put_price(S, K, maturity, r, sigma, b)
    return bond_plus_coupon - put


def capped_principal_protected_note(S, K, cap_level, maturity, r, sigma,
                                    principal, participation=1.0, n_shares=None,
                                    b=None):
    """Principal-protected note with a capped upside (a call spread).

    Like :func:`principal_protected_note` but the upside is a call spread -- long
    a call at ``K``, short a call at ``cap_level`` -- so the payoff is capped once
    the underlying passes ``cap_level``. Value is the ZC bond plus
    ``participation * n_shares * (call(K) - call(cap_level))``. At or below the
    uncapped PPN (selling the higher-strike call raises no value).
    """
    if cap_level <= K:
        raise ValueError("cap_level must exceed the strike K")
    if participation < 0:
        raise ValueError("participation must be non-negative")
    if n_shares is None:
        if S <= 0:
            raise ValueError("S must be positive to default n_shares")
        n_shares = principal / S
    bond = note_zero_coupon_bond(principal, r, maturity)
    spread = call_price(S, K, maturity, r, sigma, b) \
        - call_price(S, cap_level, maturity, r, sigma, b)
    return bond + participation * n_shares * spread


def reverse_convertible_fair_coupon(S, K, maturity, r, sigma, principal,
                                    n_shares=None, b=None, tol=1e-12,
                                    max_iter=100):
    """Coupon rate that prices a reverse convertible at par (its principal).

    Solves :func:`reverse_convertible` ``= principal`` for the ``coupon_rate``.
    The short put costs value, so the fair coupon is positive -- the enhanced yield
    that compensates the investor for the downside they sell. Closed form:
    ``coupon = (put_value / disc / principal ... )``; here solved directly since
    the note is linear in the coupon.
    """
    if n_shares is None:
        if K <= 0:
            raise ValueError("K must be positive to default n_shares")
        n_shares = principal / K
    disc = math.exp(-r * maturity)
    put = n_shares * put_price(S, K, maturity, r, sigma, b)
    # reverse_convertible = principal*(1 + c*T)*disc - put = principal
    #  => c = (principal + put - principal*disc) / (principal * disc * T)
    if maturity <= 0:
        raise ValueError("maturity must be positive")
    return (principal + put - principal * disc) / (principal * disc * maturity)


def buffered_note(S, K, buffer, maturity, r, sigma, principal, n_shares=None,
                  b=None):
    """Buffered note: absorbs the first ``buffer`` fraction of downside losses.

    The investor is short a put struck at the buffered level ``K * (1 - buffer)``
    rather than at ``K``, so losses only bite once the underlying falls more than
    ``buffer``. Value is ``PV(principal) - n_shares * put(K*(1-buffer))``. A larger
    buffer moves the put further out of the money, raising the note's value.
    """
    if not (0.0 <= buffer < 1.0):
        raise ValueError("buffer must be in [0, 1)")
    if n_shares is None:
        if K <= 0:
            raise ValueError("K must be positive to default n_shares")
        n_shares = principal / K
    buffered_strike = K * (1.0 - buffer)
    bond = note_zero_coupon_bond(principal, r, maturity)
    put = n_shares * put_price(S, buffered_strike, maturity, r, sigma, b)
    return bond - put


def note_embedded_option_value(note_value, principal, r, maturity):
    """Option component of a note: ``note_value - discounted principal``.

    Strips the guaranteed bond leg to isolate the value attributable to the
    embedded option position (positive for a PPN's long call, negative for a
    reverse convertible's short put).
    """
    return note_value - note_zero_coupon_bond(principal, r, maturity)


def phoenix_autocall_mc(S, t, r, sigma, observation_times, autocall_barrier,
                        coupon_barrier, coupon, protection_barrier=None,
                        principal=1.0, memory=True, q=0.0, n_paths=40000,
                        seed=1234567):
    """Monte Carlo a Phoenix autocallable note.

    A Phoenix pays a coupon at each observation where the spot is at or above the
    ``coupon_barrier`` (typically below the autocall level). With ``memory=True``
    any coupons missed while below the barrier are paid retroactively the next time
    the barrier is met (snowball/memory feature). If the spot reaches
    ``autocall_barrier`` the note redeems early at par plus the coupon due. At
    maturity, unredeemed, the holder gets par unless the spot is below
    ``protection_barrier`` (down-and-in), taking the downside ``principal * S_T/S``.

    Returns the discounted Monte Carlo price. Memory raises the value versus no
    memory; a lower coupon barrier pays more often. Pure standard library.
    """
    import math

    if S <= 0 or t <= 0 or sigma <= 0 or principal <= 0:
        raise ValueError("S, t, sigma, principal must be positive")
    if not observation_times or observation_times[-1] <= 0:
        raise ValueError("observation_times must be positive and non-empty")
    b = r - q
    obs = list(observation_times)
    m = len(obs)

    def _lcg(seed):
        state = seed & 0x7FFFFFFF
        cache = []

        def nxt():
            nonlocal state
            if cache:
                return cache.pop()
            state = (1103515245 * state + 12345) & 0x7FFFFFFF
            u1 = (state + 0.5) / 0x80000000
            state = (1103515245 * state + 12345) & 0x7FFFFFFF
            u2 = (state + 0.5) / 0x80000000
            rr = math.sqrt(-2.0 * math.log(u1))
            cache.append(rr * math.sin(2.0 * math.pi * u2))
            return rr * math.cos(2.0 * math.pi * u2)
        return nxt

    gen = _lcg(seed)
    total = 0.0
    for _ in range(n_paths):
        spot = S
        prev_t = 0.0
        pv = 0.0
        unpaid = 0                       # missed coupons held in memory
        redeemed = False
        for k in range(m):
            dt = obs[k] - prev_t
            prev_t = obs[k]
            z = gen()
            spot *= math.exp((b - 0.5 * sigma * sigma) * dt + sigma * math.sqrt(dt) * z)
            disc = math.exp(-r * obs[k])
            if spot >= coupon_barrier:
                due = coupon * (unpaid + 1 if memory else 1)
                pv += principal * due * disc
                unpaid = 0
            elif memory:
                unpaid += 1
            if spot >= autocall_barrier:
                pv += principal * disc     # par redemption
                redeemed = True
                break
        if not redeemed:
            disc = math.exp(-r * obs[-1])
            if protection_barrier is not None and spot < protection_barrier:
                pv += principal * (spot / S) * disc
            else:
                pv += principal * disc
        total += pv
    return total / n_paths
