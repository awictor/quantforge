"""Convert between Black (lognormal) and Bachelier (normal) implied volatilities.

Rates desks quote both a lognormal (Black-76) vol and a normal (Bachelier) vol for
the same option. The two are linked by requiring both models to return the same
price: convert by pricing under one model and inverting the other. Near the money
the leading-order relation is ``sigma_N ~ sigma_B * F`` (a lognormal move of
``sigma_B`` on a forward ``F`` is a normal move of ``sigma_B * F``). Pure standard
library (builds on the Black-76 and Bachelier pricers).
"""

from .bsm import call_price
from .bachelier import bachelier_price, bachelier_implied_vol
from .implied import implied_volatility
from .bsm import OptionType


def black_to_normal_vol(forward, strike, t, sigma_black, is_call=True):
    """Convert a Black-76 (lognormal) vol to the equivalent Bachelier (normal) vol.

    Prices the option with Black-76 at ``sigma_black`` (zero rates, forward
    measure) and inverts the Bachelier model for the normal vol that reproduces it.
    """
    if forward <= 0 or strike <= 0 or t <= 0 or sigma_black <= 0:
        raise ValueError("forward, strike, t, sigma_black must be positive")
    ot = OptionType.CALL if is_call else OptionType.PUT
    price = call_price(forward, strike, t, 0.0, sigma_black, b=0.0) if is_call \
        else _black_put(forward, strike, t, sigma_black)
    return bachelier_implied_vol(price, forward, strike, t, 0.0, option_type=ot)


def normal_to_black_vol(forward, strike, t, sigma_normal, is_call=True):
    """Convert a Bachelier (normal) vol to the equivalent Black-76 (lognormal) vol.

    Prices with Bachelier at ``sigma_normal`` and inverts Black-76. Requires
    positive forward and strike (a lognormal vol is undefined otherwise).
    """
    if forward <= 0 or strike <= 0 or t <= 0 or sigma_normal <= 0:
        raise ValueError("forward, strike, t, sigma_normal must be positive")
    ot = OptionType.CALL if is_call else OptionType.PUT
    price = bachelier_price(forward, strike, t, 0.0, sigma_normal, option_type=ot)
    # Invert Black-76 (forward measure: r=0, b=0).
    return implied_volatility(price, forward, strike, t, 0.0, b=0.0,
                              option_type=ot)


def _black_put(F, K, t, sigma):
    c = call_price(F, K, t, 0.0, sigma, b=0.0)
    return c - (F - K)          # forward-measure put-call parity, r=0
