"""Bond futures: conversion factors, delivery basis, and cheapest-to-deliver.

A bond futures contract lets the short deliver any bond from a basket, each scaled
by a *conversion factor* that normalizes it to the notional coupon. The invoice
the short receives is ``futures_price * conversion_factor + accrued``. The
delivery is optimized by choosing the *cheapest-to-deliver* bond -- the one with
the smallest net basis. This module computes the conversion factor, the invoice
price, gross and net basis, the implied repo rate, and the CTD selection. Pure
standard library.
"""

import math


def conversion_factor(coupon_rate, years_to_maturity, notional_coupon=0.06,
                      freq=2):
    """Conversion factor: price per unit face of a bond yielding the notional coupon.

    The bond's price discounted at the exchange's notional coupon
    (``notional_coupon``, e.g. 6% for CBOT Treasury futures), per unit face:

        CF = sum_i (c/freq) / (1 + y/freq)^i + 1 / (1 + y/freq)^n,

    with ``c = coupon_rate``, ``y = notional_coupon``, ``n = years * freq``. Above
    one for bonds with coupons over the notional, below one otherwise; exactly one
    when the coupon equals the notional.
    """
    if coupon_rate < 0 or years_to_maturity <= 0 or freq < 1:
        raise ValueError("invalid coupon, maturity, or frequency")
    n = int(round(years_to_maturity * freq))
    y = notional_coupon / freq
    c = coupon_rate / freq
    price = sum(c / (1.0 + y) ** i for i in range(1, n + 1)) + 1.0 / (1.0 + y) ** n
    return price


def invoice_price(futures_price, conversion_factor_, accrued_interest=0.0):
    """Invoice price the short receives: ``futures * CF + accrued``."""
    return futures_price * conversion_factor_ + accrued_interest


def gross_basis(bond_price, futures_price, conversion_factor_):
    """Gross basis ``bond_price - futures_price * CF``.

    The raw richness of the cash bond over its futures-implied forward price
    (before carry). Positive means the bond trades above its converted futures
    price.
    """
    return bond_price - futures_price * conversion_factor_


def net_basis(bond_price, futures_price, conversion_factor_, carry):
    """Net basis: gross basis less the carry to delivery.

    ``gross_basis - carry`` where ``carry`` is coupon income minus financing over
    the holding period. The cheapest-to-deliver bond minimizes the net basis.
    """
    return gross_basis(bond_price, futures_price, conversion_factor_) - carry


def implied_repo_rate(bond_price, accrued_now, futures_price,
                      conversion_factor_, accrued_delivery, days, day_count=360):
    """Implied repo rate of a deliverable bond.

    The financing rate that makes buying the bond and delivering into the future
    break even:

        repo = (invoice - dirty_now) / dirty_now * (day_count / days),

    with ``invoice = futures * CF + accrued_delivery`` and ``dirty_now =
    bond_price + accrued_now``. The cheapest-to-deliver bond has the highest
    implied repo.
    """
    if days <= 0:
        raise ValueError("days must be positive")
    dirty_now = bond_price + accrued_now
    if dirty_now <= 0:
        raise ValueError("dirty price must be positive")
    invoice = futures_price * conversion_factor_ + accrued_delivery
    return (invoice - dirty_now) / dirty_now * (day_count / days)


def cheapest_to_deliver(bonds, futures_price):
    """Index of the cheapest-to-deliver bond (minimum net basis).

    ``bonds`` is a list of dicts with keys ``price``, ``cf`` (conversion factor),
    and ``carry``. Returns the index minimizing :func:`net_basis`. Ties break to
    the first.
    """
    if not bonds:
        raise ValueError("need at least one deliverable bond")
    best_i = 0
    best_nb = None
    for i, b in enumerate(bonds):
        nb = net_basis(b["price"], futures_price, b["cf"], b.get("carry", 0.0))
        if best_nb is None or nb < best_nb:
            best_nb = nb
            best_i = i
    return best_i
