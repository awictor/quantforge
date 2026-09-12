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


def note_embedded_option_value(note_value, principal, r, maturity):
    """Option component of a note: ``note_value - discounted principal``.

    Strips the guaranteed bond leg to isolate the value attributable to the
    embedded option position (positive for a PPN's long call, negative for a
    reverse convertible's short put).
    """
    return note_value - note_zero_coupon_bond(principal, r, maturity)
