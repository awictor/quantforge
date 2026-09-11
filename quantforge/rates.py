"""Interest-rate caps, floors, and collars priced on the Bachelier (normal) model.

A **caplet** pays ``tau * max(F - K, 0)`` at the end of an accrual period, where
``F`` is the forward rate over that period, ``K`` the cap strike, and ``tau`` the
year fraction (accrual). It is a call on the forward rate, discounted by the
period's zero-coupon bond. A **cap** is a strip of caplets across successive
periods; a **floor** is the strip of floorlets (puts); a **collar** is long a
cap and short a floor.

Rates can be negative, so the caplets use the Bachelier (normal-vol) model
rather than lognormal Black-76. Each period is described by its forward rate,
its normal volatility, its accrual, its expiry (when the rate sets), and its
discount factor to the payment date.
"""

import math
from dataclasses import dataclass
from typing import Sequence, List

from .bachelier import (
    bachelier_price, bachelier_delta, bachelier_gamma, bachelier_vega,
    bachelier_implied_vol,
)
from .bsm import (
    OptionType, price as bsm_price, delta as bsm_delta, gamma as bsm_gamma,
    vega as bsm_vega,
)


@dataclass(frozen=True)
class CapletPeriod:
    forward: float          # forward rate for the period
    expiry: float           # time (years) until the rate is observed/set
    accrual: float          # year fraction tau of the accrual period
    discount: float         # discount factor to the payment date
    sigma_n: float          # normal (absolute) volatility of the forward rate


def compounded_overnight_rate(fixings, accruals) -> float:
    """Annualized rate from daily-compounding overnight fixings (SOFR-style).

    The compounded setting-in-arrears rate over a period: multiply the daily
    growth factors ``(1 + r_i tau_i)`` and annualize by the total accrual,

        rate = (prod_i (1 + r_i tau_i) - 1) / sum_i tau_i.

    ``fixings`` are the per-day annualized overnight rates and ``accruals`` the
    day-count fractions (typically 1/360). This is how compounded SOFR / SONIA
    coupons are computed.
    """
    if len(fixings) != len(accruals):
        raise ValueError("fixings and accruals must have equal length")
    if not fixings:
        raise ValueError("need at least one fixing")
    growth = 1.0
    total = 0.0
    for r, tau in zip(fixings, accruals):
        growth *= (1.0 + r * tau)
        total += tau
    if total <= 0.0:
        raise ValueError("total accrual must be positive")
    return (growth - 1.0) / total


def compounded_rate_with_lookback(fixings, accruals, lookback=0) -> float:
    """Compounded overnight rate with a lookback (observation-shift) of ``k`` days.

    Each accrual period uses the fixing observed ``lookback`` business days
    earlier (an *observation shift* also shifts the weighting to the earlier
    day). This gives the payment-lag convention used to publish a compounded SOFR
    coupon a few days before period end. ``lookback = 0`` reduces to
    :func:`compounded_overnight_rate`.

    ``fixings`` must extend at least ``lookback`` days before the accrual start;
    fixing ``i`` here is the rate applied to accrual ``i``, already shifted by the
    caller when constructing the arrays -- so this compounds ``fixings[i]`` over
    ``accruals[i]`` with an index offset used only to validate coverage.
    """
    if lookback < 0:
        raise ValueError("lookback must be non-negative")
    if len(fixings) != len(accruals):
        raise ValueError("fixings and accruals must have equal length")
    if len(fixings) <= lookback:
        raise ValueError("need more fixings than the lookback")
    growth = 1.0
    total = 0.0
    for i in range(len(accruals)):
        r = fixings[max(0, i - lookback)]
        growth *= (1.0 + r * accruals[i])
        total += accruals[i]
    if total <= 0.0:
        raise ValueError("total accrual must be positive")
    return (growth - 1.0) / total


def compounded_rate_with_lockout(fixings, accruals, lockout=0) -> float:
    """Compounded overnight rate with a rate lockout of ``k`` days.

    The final ``lockout`` business days of the period reuse the last observed
    fixing (the rate is *locked* before period end so the coupon is known early),
    the convention used for compounded fed funds. ``lockout = 0`` reduces to
    :func:`compounded_overnight_rate`.
    """
    if lockout < 0:
        raise ValueError("lockout must be non-negative")
    if len(fixings) != len(accruals):
        raise ValueError("fixings and accruals must have equal length")
    n = len(fixings)
    if n <= lockout:
        raise ValueError("need more fixings than the lockout")
    frozen = fixings[n - lockout - 1]
    growth = 1.0
    total = 0.0
    for i in range(n):
        r = fixings[i] if i < n - lockout else frozen
        growth *= (1.0 + r * accruals[i])
        total += accruals[i]
    if total <= 0.0:
        raise ValueError("total accrual must be positive")
    return (growth - 1.0) / total


def simple_average_rate(fixings, accruals) -> float:
    """Accrual-weighted arithmetic average of overnight fixings (Fed-funds style).

    ``sum_i r_i tau_i / sum_i tau_i`` -- the simple (non-compounded) averaging
    convention. Lies below the compounded rate when fixings are positive
    (compounding adds interest-on-interest).
    """
    if len(fixings) != len(accruals):
        raise ValueError("fixings and accruals must have equal length")
    if not fixings:
        raise ValueError("need at least one fixing")
    total = sum(accruals)
    if total <= 0.0:
        raise ValueError("total accrual must be positive")
    return sum(r * tau for r, tau in zip(fixings, accruals)) / total


def caplet_price(period: CapletPeriod, strike: float, is_cap: bool = True) -> float:
    """Price a single caplet (cap) or floorlet (floor).

    Value = discount * accrual * Bachelier(F, K, expiry, r=0, sigma_n),
    with the option being a call for a caplet and a put for a floorlet. The
    Bachelier price is taken undiscounted (r=0) and discounted explicitly by the
    period's bond factor, which is the market convention.
    """
    ot = OptionType.CALL if is_cap else OptionType.PUT
    if period.expiry <= 0:
        intrinsic = (max(period.forward - strike, 0.0) if is_cap
                     else max(strike - period.forward, 0.0))
        return period.discount * period.accrual * intrinsic
    undiscounted = bachelier_price(period.forward, strike, period.expiry, 0.0,
                                   period.sigma_n, ot)
    return period.discount * period.accrual * undiscounted


def caplet_greeks(period: CapletPeriod, strike: float, is_cap: bool = True):
    """Analytic Greeks of a single caplet/floorlet (normal model).

    The value is ``discount * accrual * Bachelier(F, K, expiry, 0, sigma_n)``, so
    its rate Greeks are the Bachelier Greeks in the forward rate scaled by the
    same ``discount * accrual`` factor: ``rate_delta`` (dV/dF), ``rate_gamma``
    (d2V/dF2), and ``vega`` (dV/dsigma_n). A caplet is a call on the forward, so
    its rate delta is positive; a floorlet's is negative. Returns a dict with
    ``price``, ``rate_delta``, ``rate_gamma``, ``vega``.
    """
    ot = OptionType.CALL if is_cap else OptionType.PUT
    f = period.forward
    scale = period.discount * period.accrual
    price = caplet_price(period, strike, is_cap)
    if period.expiry <= 0:
        itm = (f > strike) if is_cap else (f < strike)
        rate_delta = scale * ((1.0 if is_cap else -1.0) if itm else 0.0)
        return {"price": price, "rate_delta": rate_delta, "rate_gamma": 0.0,
                "vega": 0.0}
    rate_delta = scale * bachelier_delta(f, strike, period.expiry, 0.0,
                                         period.sigma_n, ot)
    rate_gamma = scale * bachelier_gamma(f, strike, period.expiry, 0.0,
                                         period.sigma_n)
    vega = scale * bachelier_vega(f, strike, period.expiry, 0.0, period.sigma_n)
    return {"price": price, "rate_delta": rate_delta, "rate_gamma": rate_gamma,
            "vega": vega}


def _sum_greeks(periods, strike, is_cap):
    out = {"price": 0.0, "rate_delta": 0.0, "rate_gamma": 0.0, "vega": 0.0}
    for p in periods:
        g = caplet_greeks(p, strike, is_cap)
        for k in out:
            out[k] += g[k]
    return out


def cap_price(periods: Sequence[CapletPeriod], strike: float) -> float:
    """Price an interest-rate cap as the sum of its caplets."""
    return sum(caplet_price(p, strike, is_cap=True) for p in periods)


def cap_greeks(periods: Sequence[CapletPeriod], strike: float):
    """Aggregate Greeks of a cap: the summed caplet ``price``/``rate_delta``/
    ``rate_gamma``/``vega`` (all per unit notional)."""
    return _sum_greeks(periods, strike, True)


def floor_greeks(periods: Sequence[CapletPeriod], strike: float):
    """Aggregate Greeks of a floor: the summed floorlet Greeks."""
    return _sum_greeks(periods, strike, False)


def floor_price(periods: Sequence[CapletPeriod], strike: float) -> float:
    """Price an interest-rate floor as the sum of its floorlets."""
    return sum(caplet_price(p, strike, is_cap=False) for p in periods)


def _sabr_caplet_price(period, strike, alpha, beta, rho, nu, is_cap):
    """One caplet/floorlet priced at its own SABR normal vol.

    Each caplet sets on its forward ``period.forward`` at ``period.expiry``, so
    the SABR normal vol is read at that ``(F, K, expiry)`` point and used in the
    Bachelier caplet formula (rates can be negative, hence the normal model).
    """
    from .sabr import sabr_normal_vol
    ot = OptionType.CALL if is_cap else OptionType.PUT
    if period.expiry <= 0:
        intrinsic = (max(period.forward - strike, 0.0) if is_cap
                     else max(strike - period.forward, 0.0))
        return period.discount * period.accrual * intrinsic
    vol = sabr_normal_vol(period.forward, strike, period.expiry,
                          alpha, beta, rho, nu)
    undiscounted = bachelier_price(period.forward, strike, period.expiry, 0.0,
                                   vol, ot)
    return period.discount * period.accrual * undiscounted


def caplet_implied_normal_vol(price, period: CapletPeriod, strike,
                              is_cap=True) -> float:
    """Normal (Bachelier) implied vol of a caplet/floorlet from its price.

    Divides out the ``discount * accrual`` factor to recover the undiscounted
    Bachelier option value, then inverts it with :func:`bachelier_implied_vol`.
    Inverse of :func:`caplet_price` in the forward rate.
    """
    if period.expiry <= 0:
        raise ValueError("cannot imply vol at or past expiry")
    scale = period.discount * period.accrual
    if scale <= 0:
        raise ValueError("discount * accrual must be positive")
    ot = OptionType.CALL if is_cap else OptionType.PUT
    undiscounted = price / scale
    return bachelier_implied_vol(undiscounted, period.forward, strike,
                                 period.expiry, 0.0, ot)


def sabr_cap_price(periods: Sequence[CapletPeriod], strike,
                   alpha, beta, rho, nu) -> float:
    """Price a cap under a single SABR smile (normal model).

    Sums caplets, each valued at the SABR normal vol read at its own forward and
    expiry -- so one calibrated ``(alpha, beta, rho, nu)`` prices the whole cap
    consistently across the smile, rather than a flat per-period ``sigma_n``.
    """
    return sum(_sabr_caplet_price(p, strike, alpha, beta, rho, nu, True)
               for p in periods)


def sabr_floor_price(periods: Sequence[CapletPeriod], strike,
                     alpha, beta, rho, nu) -> float:
    """Price a floor under a single SABR smile (normal model): the sum of
    floorlets, each at the SABR normal vol of its own forward and expiry."""
    return sum(_sabr_caplet_price(p, strike, alpha, beta, rho, nu, False)
               for p in periods)


def collar_price(periods: Sequence[CapletPeriod], cap_strike: float,
                 floor_strike: float) -> float:
    """Price a collar: long a cap at ``cap_strike``, short a floor at ``floor_strike``.

    The net value is ``cap - floor``; a zero-cost collar is the pair of strikes
    that makes this zero.
    """
    return cap_price(periods, cap_strike) - floor_price(periods, floor_strike)


def caplet_floorlet_parity(period: CapletPeriod, strike: float) -> float:
    """Caplet - floorlet at the same strike = discounted forward-minus-strike.

    A put-call-parity identity used to check the pricer:
    ``caplet - floorlet = discount * accrual * (F - K)``.
    """
    return period.discount * period.accrual * (period.forward - strike)


def annuity(periods: Sequence[CapletPeriod]) -> float:
    """Present-value annuity (level / PV01) of a swap: sum of accrual*discount."""
    return sum(p.accrual * p.discount for p in periods)


def swaption_price(swap_rate, strike, expiry, sigma_n, periods,
                   payer=True) -> float:
    """Bachelier price of a European swaption on the underlying swap.

    A payer swaption is a call on the swap rate; a receiver is a put. The value
    is the swap's PV annuity times a Bachelier option on the forward swap rate:

        V = annuity * Bachelier(swap_rate, strike, expiry, r=0, sigma_n).

    Args:
        swap_rate: current forward swap rate.
        strike: fixed strike rate.
        expiry: option expiry (years) — when the swap rate sets.
        sigma_n: normal (absolute) volatility of the swap rate.
        periods: the underlying swap's ``CapletPeriod`` legs, used only for the
            annuity (accrual and discount factors).
        payer: True for a payer (call), False for a receiver (put).

    Rates may be negative; the normal model handles that.
    """
    ann = annuity(periods)
    ot = OptionType.CALL if payer else OptionType.PUT
    if expiry <= 0:
        intrinsic = (max(swap_rate - strike, 0.0) if payer
                     else max(strike - swap_rate, 0.0))
        return ann * intrinsic
    undiscounted = bachelier_price(swap_rate, strike, expiry, 0.0, sigma_n, ot)
    return ann * undiscounted


def swaption_greeks(swap_rate, strike, expiry, sigma_n, periods, payer=True):
    """Analytic Greeks of a European swaption (normal model).

    The value is ``annuity * Bachelier(swap_rate, strike, expiry, 0, sigma_n)``,
    so its swap-rate Greeks are the Bachelier Greeks scaled by the annuity:
    ``rate_delta`` (dV/d swap_rate), ``rate_gamma`` (d2V/d swap_rate^2), and
    ``vega`` (dV/dsigma_n). A payer swaption is a call on the swap rate (positive
    rate delta); a receiver is a put (negative). Returns a dict with ``price``,
    ``rate_delta``, ``rate_gamma``, ``vega``, ``annuity``.
    """
    ann = annuity(periods)
    ot = OptionType.CALL if payer else OptionType.PUT
    price = swaption_price(swap_rate, strike, expiry, sigma_n, periods, payer)
    if expiry <= 0:
        itm = (swap_rate > strike) if payer else (swap_rate < strike)
        rate_delta = ann * ((1.0 if payer else -1.0) if itm else 0.0)
        return {"price": price, "rate_delta": rate_delta, "rate_gamma": 0.0,
                "vega": 0.0, "annuity": ann}
    rate_delta = ann * bachelier_delta(swap_rate, strike, expiry, 0.0, sigma_n, ot)
    rate_gamma = ann * bachelier_gamma(swap_rate, strike, expiry, 0.0, sigma_n)
    vega = ann * bachelier_vega(swap_rate, strike, expiry, 0.0, sigma_n)
    return {"price": price, "rate_delta": rate_delta, "rate_gamma": rate_gamma,
            "vega": vega, "annuity": ann}


def swaption_parity(swap_rate, strike, periods) -> float:
    """Payer - receiver at the same strike = annuity * (swap_rate - strike)."""
    return annuity(periods) * (swap_rate - strike)


def black_swaption_price(swap_rate, strike, expiry, sigma_b, periods,
                         payer=True) -> float:
    """Black (lognormal) price of a European swaption on the underlying swap.

    The market-standard lognormal counterpart to :func:`swaption_price`: the
    forward swap rate is modelled as lognormal with (Black) volatility
    ``sigma_b``, and the swaption is the annuity times a zero-carry Black-76
    option on the rate:

        V = annuity * Black76(swap_rate, strike, expiry, sigma_b).

    A payer swaption is a call on the rate, a receiver a put. Requires positive
    ``swap_rate`` and ``strike`` (use :func:`swaption_price` for the normal model
    when rates may be negative).
    """
    if swap_rate <= 0 or strike <= 0:
        raise ValueError("Black swaption needs positive swap_rate and strike")
    ann = annuity(periods)
    ot = OptionType.CALL if payer else OptionType.PUT
    if expiry <= 0:
        intrinsic = (max(swap_rate - strike, 0.0) if payer
                     else max(strike - swap_rate, 0.0))
        return ann * intrinsic
    # Black-76: BSM with zero cost of carry (b = 0), undiscounted (r = 0).
    undiscounted = bsm_price(swap_rate, strike, expiry, 0.0, sigma_b, ot, b=0.0)
    return ann * undiscounted


def black_swaption_greeks(swap_rate, strike, expiry, sigma_b, periods,
                          payer=True):
    """Analytic Greeks of a Black (lognormal) European swaption.

    The value is ``annuity * Black76(swap_rate, strike, expiry, sigma_b)``, so
    the swap-rate Greeks are the Black-76 Greeks scaled by the annuity:
    ``rate_delta`` (dV/d swap_rate), ``rate_gamma`` (d2V/d swap_rate^2), and
    ``vega`` (dV/dsigma_b). Returns a dict with ``price``, ``rate_delta``,
    ``rate_gamma``, ``vega``, ``annuity``.
    """
    if swap_rate <= 0 or strike <= 0:
        raise ValueError("Black swaption needs positive swap_rate and strike")
    ann = annuity(periods)
    ot = OptionType.CALL if payer else OptionType.PUT
    price = black_swaption_price(swap_rate, strike, expiry, sigma_b, periods, payer)
    if expiry <= 0:
        itm = (swap_rate > strike) if payer else (swap_rate < strike)
        rate_delta = ann * ((1.0 if payer else -1.0) if itm else 0.0)
        return {"price": price, "rate_delta": rate_delta, "rate_gamma": 0.0,
                "vega": 0.0, "annuity": ann}
    rate_delta = ann * bsm_delta(swap_rate, strike, expiry, 0.0, sigma_b, ot, b=0.0)
    rate_gamma = ann * bsm_gamma(swap_rate, strike, expiry, 0.0, sigma_b, b=0.0)
    vega = ann * bsm_vega(swap_rate, strike, expiry, 0.0, sigma_b, b=0.0)
    return {"price": price, "rate_delta": rate_delta, "rate_gamma": rate_gamma,
            "vega": vega, "annuity": ann}


def sabr_swaption_price(swap_rate, strike, expiry, periods,
                        alpha, beta, rho, nu, payer=True, model="black"):
    """Price a European swaption whose smile is described by a SABR model.

    Reads the SABR-implied volatility at the (forward swap rate, strike, expiry)
    point and feeds it into the matching swaption pricer:

      * ``model="black"``  -> Hagan lognormal vol :func:`quantforge.sabr_vol`
        into :func:`black_swaption_price` (requires positive rate and strike);
      * ``model="normal"`` -> Hagan normal vol :func:`quantforge.sabr_normal_vol`
        into :func:`swaption_price` (handles negative rates).

    This is the standard way SABR is used on swaptions: one calibrated smile
    prices every strike consistently. Returns the swaption present value.
    """
    from .sabr import sabr_vol, sabr_normal_vol
    if model == "black":
        vol = sabr_vol(swap_rate, strike, expiry, alpha, beta, rho, nu)
        return black_swaption_price(swap_rate, strike, expiry, vol, periods, payer)
    if model == "normal":
        vol = sabr_normal_vol(swap_rate, strike, expiry, alpha, beta, rho, nu)
        return swaption_price(swap_rate, strike, expiry, vol, periods, payer)
    raise ValueError("model must be 'black' or 'normal'")


def swaption_implied_normal_vol(price, swap_rate, strike, expiry, periods,
                                payer=True) -> float:
    """Normal (Bachelier) implied vol of a swaption from its price.

    Divides out the annuity to recover the undiscounted Bachelier option value,
    then inverts it with :func:`bachelier_implied_vol`. Inverse of
    :func:`swaption_price`.
    """
    ann = annuity(periods)
    if ann <= 0:
        raise ValueError("annuity must be positive")
    if expiry <= 0:
        raise ValueError("cannot imply vol at or past expiry")
    ot = OptionType.CALL if payer else OptionType.PUT
    undiscounted = price / ann
    return bachelier_implied_vol(undiscounted, swap_rate, strike, expiry, 0.0, ot)


def swaption_implied_black_vol(price, swap_rate, strike, expiry, periods,
                               payer=True) -> float:
    """Black (lognormal) implied vol of a swaption from its price.

    Divides out the annuity to recover the undiscounted Black-76 option value,
    then inverts it with :func:`quantforge.implied_volatility` at zero carry.
    Inverse of :func:`black_swaption_price`; requires positive rate and strike.
    """
    from .implied import implied_volatility
    if swap_rate <= 0 or strike <= 0:
        raise ValueError("Black swaption needs positive swap_rate and strike")
    ann = annuity(periods)
    if ann <= 0:
        raise ValueError("annuity must be positive")
    if expiry <= 0:
        raise ValueError("cannot imply vol at or past expiry")
    ot = OptionType.CALL if payer else OptionType.PUT
    undiscounted = price / ann
    return implied_volatility(undiscounted, swap_rate, strike, expiry, 0.0, ot,
                              b=0.0)
