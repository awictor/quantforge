"""Warrants and employee stock options (ESOs) with dilution and forfeiture.

Warrants issued by the company create new shares on exercise, diluting existing
holders, so their value is the vanilla call scaled by ``M / (M + N)`` for ``M``
existing and ``N`` new (warrant) shares. Employee stock options are long-dated
calls whose holders exercise early and forfeit before vesting; the FASB-endorsed
practical model (Hull-White) values them at an expected life shorter than the
contractual term and multiplies by the probability of surviving to exercise given
a post-vest exit rate. Pure standard library on top of :mod:`quantforge.bsm`.
"""

import math

from .bsm import call_price


def conversion_value(S, conversion_ratio):
    """Parity (conversion) value of a convertible: ``conversion_ratio * S``.

    The worth of the shares the bond converts into -- the equity floor of the
    convertible.
    """
    if conversion_ratio < 0:
        raise ValueError("conversion_ratio must be non-negative")
    return conversion_ratio * S


def straight_bond_floor(face, coupon_rate, maturity, r, credit_spread=0.0,
                        freq=2):
    """Investment (bond) value of a convertible ignoring the conversion option.

    Discounts the straight bond's coupons and principal at the risk-free rate plus
    a ``credit_spread`` (continuously compounded). This is the debt floor: the
    convertible cannot be worth less than this if held to maturity without
    converting.
    """
    if face <= 0 or maturity <= 0 or freq < 1:
        raise ValueError("face, maturity must be positive and freq >= 1")
    y = r + credit_spread
    n = int(round(maturity * freq))
    cpn = face * coupon_rate / freq
    pv = 0.0
    for i in range(1, n + 1):
        t = i / freq
        cash = cpn + (face if i == n else 0.0)
        pv += cash * math.exp(-y * t)
    return pv


def convertible_bond_value(S, conversion_ratio, face, coupon_rate, maturity, r,
                           sigma, credit_spread=0.0, freq=2):
    """Convertible bond value via the component (bond floor + call) approximation.

    Values the convertible as its :func:`straight_bond_floor` plus a call option on
    the ``conversion_ratio`` shares struck at the floor's per-share equivalent --
    the standard decomposition ``CB = bond floor + conversion_ratio *
    call(S, K=face/ratio adjusted)``. Here the call strike is set so that at
    maturity the holder converts when ``conversion_ratio * S > face``, i.e. strike
    ``= face / conversion_ratio`` on ``conversion_ratio`` shares. The result is at
    least the bond floor and at least the conversion value.
    """
    if conversion_ratio <= 0:
        raise ValueError("conversion_ratio must be positive")
    floor = straight_bond_floor(face, coupon_rate, maturity, r, credit_spread, freq)
    strike = face / conversion_ratio
    option = conversion_ratio * call_price(S, strike, maturity, r, sigma, b=r)
    return floor + option


def conversion_premium(convertible_price, S, conversion_ratio):
    """Conversion (equity) premium: how far the convertible trades above parity.

    ``convertible_price / conversion_value - 1`` -- the fractional premium an
    investor pays over the value of the underlying shares for the bond's downside
    protection. Non-negative when the convertible trades at or above parity.
    """
    parity = conversion_value(S, conversion_ratio)
    if parity <= 0:
        raise ValueError("conversion value must be positive")
    return convertible_price / parity - 1.0


def investment_premium(convertible_price, bond_floor):
    """Investment premium: how far the convertible trades above its bond floor.

    ``convertible_price / bond_floor - 1`` -- the fractional premium over the
    straight-debt value, paid for the equity upside. Non-negative when the
    convertible trades at or above its floor.
    """
    if bond_floor <= 0:
        raise ValueError("bond_floor must be positive")
    return convertible_price / bond_floor - 1.0


def convertible_breakeven_years(convertible_price, S, conversion_ratio,
                                bond_coupon_income, dividend_income):
    """Years for extra income to recoup the conversion premium.

    ``(convertible_price - conversion_value) / (bond_coupon_income -
    dividend_income)`` -- the time for the convertible's income advantage over the
    equivalent shares to pay back the dollar conversion premium. Requires the bond
    to yield more than the shares (positive net income); returns ``inf`` if not.
    """
    premium_dollars = convertible_price - conversion_value(S, conversion_ratio)
    net_income = bond_coupon_income - dividend_income
    if net_income <= 0.0:
        return float("inf")
    return premium_dollars / net_income


def pv_dividends(dividends, r):
    """Present value of a discrete dividend schedule ``[(t, amount), ...]``.

    Each cash dividend is discounted at the continuously-compounded rate ``r``:
    ``sum_i D_i e^{-r t_i}``. Dividends at or before time zero (``t <= 0``) are
    excluded (already paid).
    """
    total = 0.0
    for t, d in dividends:
        if d < 0:
            raise ValueError("dividend amounts must be non-negative")
        if t > 0:
            total += d * math.exp(-r * t)
    return total


def discrete_dividend_price(S, K, t, r, sigma, dividends, is_call=True):
    """European option on a stock paying known discrete cash dividends.

    Uses the escrowed-dividend (spot-minus-PV-of-dividends) approximation: the
    option is priced with Black-Scholes on the dividend-adjusted spot
    ``S - PV(dividends up to expiry)`` and carry ``b = r`` (the adjusted spot grows
    at the risk-free rate). Reduces to the plain BSM call/put when there are no
    dividends before expiry. Put and call satisfy
    ``C - P = (S - PV_div) - K e^{-r t}``.
    """
    from .bsm import put_price
    pv = pv_dividends([(dt, d) for dt, d in dividends if dt <= t], r)
    s_adj = S - pv
    if s_adj <= 0:
        raise ValueError("dividend PV exceeds spot; adjusted spot non-positive")
    if is_call:
        return call_price(s_adj, K, t, r, sigma, b=r)
    return put_price(s_adj, K, t, r, sigma, b=r)


def forward_with_dividends(S, t, r, dividends):
    """Forward price of a stock paying discrete dividends: ``(S - PV_div) e^{r t}``.

    The escrowed-dividend forward: the dividend-stripped spot compounded at the
    risk-free rate. Equals ``S e^{r t}`` when no dividends fall before ``t``.
    """
    pv = pv_dividends([(dt, d) for dt, d in dividends if dt <= t], r)
    return (S - pv) * math.exp(r * t)


def dilution_factor(existing_shares, new_shares):
    """Dilution multiplier ``M / (M + N)`` for issuing ``N`` new shares on ``M``.

    Below one whenever new shares are created; one when ``N = 0``. Applied to a
    warrant's per-option payoff because exercise expands the share count.
    """
    if existing_shares <= 0:
        raise ValueError("existing_shares must be positive")
    if new_shares < 0:
        raise ValueError("new_shares must be non-negative")
    return existing_shares / (existing_shares + new_shares)


def warrant_price(S, K, t, r, sigma, existing_shares, new_shares, b=None):
    """Warrant value: the vanilla call scaled by the :func:`dilution_factor`.

    A single-period dilution adjustment -- the standard textbook approximation
    ``value = M/(M+N) * call(S, K, ...)`` -- which reduces to the plain call when
    no new shares are issued. ``S`` is the current (pre-dilution) share price.
    """
    factor = dilution_factor(existing_shares, new_shares)
    return factor * call_price(S, K, t, r, sigma, b)


def eso_expected_life(vesting, contractual_term, exit_rate):
    """Expected life of an ESO given post-vest exit and the contractual term.

    After vesting, holders leave (and exercise or forfeit) at a constant hazard
    ``exit_rate``; the expected time to exercise, capped at the contractual term,
    is

        vesting + (1 - e^{-exit_rate * (T - vesting)}) / exit_rate

    -- between ``vesting`` (immediate exit) and ``contractual_term`` (no exit).
    """
    if not (0.0 <= vesting <= contractual_term):
        raise ValueError("vesting must lie in [0, contractual_term]")
    if exit_rate < 0:
        raise ValueError("exit_rate must be non-negative")
    post = contractual_term - vesting
    if exit_rate == 0.0:
        return contractual_term
    return vesting + (1.0 - math.exp(-exit_rate * post)) / exit_rate


def eso_value(S, K, contractual_term, r, sigma, vesting, exit_rate,
              forfeiture_rate=0.0, b=None):
    """Employee stock option value (Hull-White practical / FASB 123R style).

    Prices the ESO as a call at the :func:`eso_expected_life` rather than the full
    term, then multiplies by the probability of surviving pre-vest forfeiture
    ``e^{-forfeiture_rate * vesting}``. Cheaper than the vanilla call on the
    contractual term because early exercise shortens the option and forfeiture
    can extinguish it. With ``exit_rate = 0`` and ``forfeiture_rate = 0`` it
    reduces to the vanilla call.
    """
    if forfeiture_rate < 0:
        raise ValueError("forfeiture_rate must be non-negative")
    life = eso_expected_life(vesting, contractual_term, exit_rate)
    survive = math.exp(-forfeiture_rate * vesting)
    return survive * call_price(S, K, life, r, sigma, b)
