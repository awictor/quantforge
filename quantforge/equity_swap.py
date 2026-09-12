"""Equity swaps: total-return swaps and dividend swaps.

A total-return swap pays the equity total return (price change plus dividends) on
a notional against a financing leg (a funding rate plus spread on the same
notional). A dividend swap exchanges the realized dividends of an index over a
period for a fixed strike; its fair strike is the present value of the expected
dividends grossed to the settlement convention. This module values both legs and
solves the fair terms. Pure standard library.
"""

import math


def total_return_leg(notional, start_price, end_price, dividends):
    """Equity total-return leg: price return plus dividends on the notional.

    ``notional * ((end_price - start_price + dividends) / start_price)`` -- the
    cash the total-return receiver collects (negative if the equity fell more than
    its dividends).
    """
    if start_price <= 0:
        raise ValueError("start_price must be positive")
    return notional * ((end_price - start_price + dividends) / start_price)


def financing_leg(notional, funding_rate, spread, year_fraction):
    """Financing leg of a TRS: ``notional * (funding_rate + spread) * tau``.

    The interest the total-return receiver pays on the notional over the accrual
    period ``year_fraction``.
    """
    return notional * (funding_rate + spread) * year_fraction


def total_return_swap_value(notional, start_price, end_price, dividends,
                            funding_rate, spread, year_fraction):
    """Net value to the total-return receiver: equity leg minus financing leg.

    Positive when the equity total return beats the financing cost.
    """
    equity = total_return_leg(notional, start_price, end_price, dividends)
    funding = financing_leg(notional, funding_rate, spread, year_fraction)
    return equity - funding


def trs_fair_spread(start_price, expected_end_price, expected_dividends,
                    funding_rate, year_fraction):
    """Financing spread that zeroes the expected TRS value.

    Solves ``E[equity return] = (funding_rate + spread) * tau`` for the spread:

        spread = E[total return] / tau - funding_rate,

    with ``E[total return] = (E[end] - start + E[div]) / start``. The spread the
    financing leg must carry so the swap is fair at inception.
    """
    if start_price <= 0:
        raise ValueError("start_price must be positive")
    if year_fraction <= 0:
        raise ValueError("year_fraction must be positive")
    total_return = (expected_end_price - start_price + expected_dividends) / start_price
    return total_return / year_fraction - funding_rate


def dividend_swap_fair_strike(expected_dividends, discount_factors):
    """Fair strike of a dividend swap: PV of the expected dividend stream.

    ``sum_i DF_i * D_i`` -- the present value of the expected dividends the
    floating leg will pay, which the fixed strike must match for the swap to have
    zero value at inception.
    """
    if len(expected_dividends) != len(discount_factors):
        raise ValueError("expected_dividends and discount_factors must align")
    return sum(df * d for d, df in zip(expected_dividends, discount_factors))


def dividend_swap_value(realized_dividends, strike, discount_factors,
                        notional=1.0):
    """Value of a dividend swap to the fixed-strike payer (dividend receiver).

    ``notional * (sum_i DF_i * D_i - strike)`` -- the PV of realized dividends less
    the fixed strike. Zero when the strike equals the
    :func:`dividend_swap_fair_strike` for the realized stream.
    """
    if len(realized_dividends) != len(discount_factors):
        raise ValueError("realized_dividends and discount_factors must align")
    pv = sum(df * d for d, df in zip(realized_dividends, discount_factors))
    return notional * (pv - strike)
