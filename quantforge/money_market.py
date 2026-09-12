"""Money-market yield conventions: discount, CD, and bond-equivalent yields.

Short-term instruments (T-bills, commercial paper, CDs) quote on different day-
count and compounding conventions. This module converts between the bank discount
yield, the price, the money-market (CD/actual-360) yield, and the bond-equivalent
yield (actual-365, comparable to a coupon bond). Pure standard library.
"""


def price_from_discount(face, discount_rate, days, year_days=360):
    """Price of a discount instrument from its bank discount rate.

    ``price = face * (1 - discount_rate * days / year_days)`` -- the bank discount
    convention prices off the *face*, not the price, on an actual/360 basis.
    """
    if face <= 0 or days <= 0:
        raise ValueError("face and days must be positive")
    return face * (1.0 - discount_rate * days / year_days)


def bank_discount_yield(face, price, days, year_days=360):
    """Bank discount yield from the price: ``(face - price)/face * year_days/days``.

    Inverse of :func:`price_from_discount`. Understates the true return because it
    divides the discount by face rather than price.
    """
    if face <= 0 or price <= 0 or days <= 0:
        raise ValueError("face, price, days must be positive")
    return (face - price) / face * (year_days / days)


def money_market_yield(face, price, days, year_days=360):
    """Money-market (CD-equivalent) yield: return on *price*, actual/360.

    ``(face - price)/price * year_days/days`` -- the actual return per invested
    dollar, higher than the :func:`bank_discount_yield` (which divides by face).
    """
    if face <= 0 or price <= 0 or days <= 0:
        raise ValueError("face, price, days must be positive")
    return (face - price) / price * (year_days / days)


def bond_equivalent_yield(face, price, days, year_days=365):
    """Bond-equivalent yield: return on price on an actual/365 basis.

    ``(face - price)/price * 365/days`` -- puts a discount instrument on the same
    footing as a coupon bond (actual/365), so it exceeds both the discount and the
    actual/360 money-market yield.
    """
    if face <= 0 or price <= 0 or days <= 0:
        raise ValueError("face, price, days must be positive")
    return (face - price) / price * (year_days / days)


def discount_to_bond_equivalent(discount_rate, days):
    """Convert a bank discount rate directly to a bond-equivalent yield.

    Prices at par 100 off the discount, then takes the actual/365 return on price.
    Always above the input discount rate.
    """
    price = price_from_discount(100.0, discount_rate, days)
    return bond_equivalent_yield(100.0, price, days)


def holding_period_return(buy_price, sell_price, income=0.0):
    """Holding-period return ``(sell - buy + income) / buy``.

    The total return over the holding period including any interim income; not
    annualized.
    """
    if buy_price <= 0:
        raise ValueError("buy_price must be positive")
    return (sell_price - buy_price + income) / buy_price
