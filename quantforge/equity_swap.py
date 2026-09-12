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


def variance_swap_payoff(realized_vol, strike_vol, variance_notional):
    """Variance-swap payoff ``variance_notional * (realized_vol^2 - strike_vol^2)``.

    Settles on the difference between realized and strike *variance* (vols entered
    as decimals, e.g. 0.20 for 20%). Convex in realized vol -- the variance
    convention penalizes large moves more than a volatility swap.
    """
    if variance_notional < 0:
        raise ValueError("variance_notional must be non-negative")
    return variance_notional * (realized_vol ** 2 - strike_vol ** 2)


def vega_notional_to_variance_notional(vega_notional, strike_vol):
    """Convert a vega notional to the equivalent variance notional.

    ``variance_notional = vega_notional / (2 * strike_vol)`` -- the market quotes
    variance swaps in vega terms (P&L per vol point at the strike); this is the
    variance notional that reproduces that sensitivity.
    """
    if strike_vol <= 0:
        raise ValueError("strike_vol must be positive")
    return vega_notional / (2.0 * strike_vol)


def volatility_swap_payoff(realized_vol, strike_vol, vega_notional):
    """Volatility-swap payoff ``vega_notional * (realized_vol - strike_vol)``.

    Linear in realized vol (in vol points). Unlike the variance swap it has no
    convexity, so it prices below a variance swap struck at the same vol (the
    convexity value / vol-of-vol adjustment).
    """
    if vega_notional < 0:
        raise ValueError("vega_notional must be non-negative")
    return vega_notional * (realized_vol - strike_vol)


def variance_swap_mtm(accrued_variance, expected_future_variance, elapsed,
                      total_time, strike_vol, variance_notional, discount_factor):
    """Mark-to-market of a seasoned variance swap.

    Blends the realized (accrued) variance over the elapsed fraction with the
    expected future variance over the remainder to get the expected terminal
    realized variance, then discounts the variance-swap payoff:

        E[realized_var] = (elapsed * accrued + (total - elapsed) * future) / total
        MTM = DF * variance_notional * (E[realized_var] - strike_vol^2).

    At inception (``elapsed = 0``) it is the discounted expected-minus-strike
    variance; at expiry it is the fully-realized payoff.
    """
    if total_time <= 0 or elapsed < 0 or elapsed > total_time:
        raise ValueError("require 0 <= elapsed <= total_time and total_time > 0")
    exp_var = (elapsed * accrued_variance
               + (total_time - elapsed) * expected_future_variance) / total_time
    return discount_factor * variance_notional * (exp_var - strike_vol ** 2)


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
