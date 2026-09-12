"""Dual-currency deposits (DCD): yield-enhanced FX-linked deposits.

A dual-currency deposit pays an above-market coupon in exchange for the depositor
selling an FX option: at maturity the bank may repay in the alternate currency at
a pre-agreed strike if that is cheaper for it. The enhanced yield is the base
deposit rate plus the annualized option premium the depositor earns for the sold
conversion right. This module computes the enhanced yield, the maturity payoff in
base-currency terms, and the breakeven spot. Pure standard library on top of
:mod:`quantforge.bsm`.
"""

import math

from .bsm import call_price, put_price


def dcd_enhanced_yield(base_deposit_rate, option_premium_rate, tenor):
    """Enhanced annual yield of a DCD: base rate plus annualized premium.

    ``base_deposit_rate + option_premium_rate / tenor`` -- the option premium
    (as a fraction of notional, earned once) spread over the ``tenor`` in years and
    added to the plain deposit rate. Above the base deposit rate whenever the sold
    option has value.
    """
    if tenor <= 0:
        raise ValueError("tenor must be positive")
    if option_premium_rate < 0:
        raise ValueError("option_premium_rate must be non-negative")
    return base_deposit_rate + option_premium_rate / tenor


def dcd_option_premium_rate(spot, strike, tenor, r_domestic, r_foreign, sigma,
                            deposit_ccy_is_base=True):
    """Premium rate (fraction of notional) of the option embedded in a DCD.

    The depositor sells a call on the base currency (if depositing the base
    currency) or a put (if depositing the alternate), struck at the conversion
    ``strike``, priced with Garman-Kohlhagen (carry ``b = r_domestic -
    r_foreign``). Returned per unit notional. Higher vol and a nearer strike raise
    the premium, and thus the DCD's yield pickup.
    """
    if spot <= 0 or strike <= 0 or tenor <= 0:
        raise ValueError("spot, strike, tenor must be positive")
    b = r_domestic - r_foreign
    if deposit_ccy_is_base:
        return call_price(spot, strike, tenor, r_domestic, sigma, b) / spot
    return put_price(spot, strike, tenor, r_domestic, sigma, b) / spot


def dcd_maturity_payoff(notional, coupon_rate, tenor, spot_at_maturity, strike,
                        deposit_ccy_is_base=True):
    """Base-currency value of a DCD at maturity, including conversion.

    The depositor always earns ``notional * (1 + coupon_rate * tenor)`` of
    principal-plus-coupon; if the option the depositor sold finishes in the money
    the repayment is converted at ``strike`` rather than ``spot_at_maturity``,
    costing the depositor the shortfall. Returns the base-currency value received.

    For a base-currency deposit (sold call), conversion bites when
    ``spot_at_maturity > strike``: repayment happens at the worse ``strike``, so
    value = principal_plus_coupon * strike / spot when converted.
    """
    if notional <= 0 or tenor <= 0 or spot_at_maturity <= 0 or strike <= 0:
        raise ValueError("notional, tenor, spot, strike must be positive")
    gross = notional * (1.0 + coupon_rate * tenor)
    if deposit_ccy_is_base:
        # Converted (repaid in alternate ccy at strike) when spot > strike.
        if spot_at_maturity > strike:
            return gross * strike / spot_at_maturity
        return gross
    # Alternate-ccy deposit (sold put): converted when spot < strike.
    if spot_at_maturity < strike:
        return gross * strike / spot_at_maturity
    return gross


def dcd_breakeven_spot(coupon_rate, base_deposit_rate, tenor, strike,
                       deposit_ccy_is_base=True):
    """Spot at which the DCD matches a plain deposit (breakeven).

    The converted DCD value equals the plain-deposit value when the extra coupon
    exactly offsets the conversion loss. For a base deposit (converted above
    ``strike``):

        strike * (1 + coupon*tenor) / breakeven = 1 + base*tenor
        => breakeven = strike * (1 + coupon*tenor) / (1 + base*tenor).

    Above the strike (in the converted region) for a positive yield pickup.
    """
    if tenor <= 0 or strike <= 0:
        raise ValueError("tenor and strike must be positive")
    ratio = (1.0 + coupon_rate * tenor) / (1.0 + base_deposit_rate * tenor)
    if deposit_ccy_is_base:
        return strike * ratio
    return strike / ratio
