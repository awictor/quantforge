"""Variance- and volatility-swap fair strikes via static option replication.

A variance swap pays the realized variance of the underlying over its life. Its
fair strike (the delivery variance that makes the swap worth zero at inception)
equals the risk-neutral expectation of realized variance, which Carr-Madan /
Demeterfi-Derman-Kamani-Zou show is replicated by a *strip* of options:

    K_var = (2/t) [ r t - (S0 e^{r t}/K* - 1) - ln(K*/S0)
                    + e^{r t} ( put strip below K*  +  call strip above K*) ]

where the strips integrate option prices weighted by ``1/K^2`` and ``K*`` is the
forward split point. Given a discrete option chain we evaluate the strip by the
trapezoidal rule. The fair *volatility*-swap strike is approximately
``sqrt(K_var)`` minus a convexity adjustment; we return the plain
``sqrt(K_var)`` proxy and expose the variance strike as the primary result.

This is the standard desk method for marking variance swaps directly from
listed option prices, model-free.
"""

import math
from typing import Sequence, Callable


def _strip_integral(strikes, prices, weight):
    """Trapezoidal integral of ``weight(K) * price(K)`` over the given strikes."""
    total = 0.0
    for i in range(len(strikes) - 1):
        k0, k1 = strikes[i], strikes[i + 1]
        g0 = weight(k0) * prices[i]
        g1 = weight(k1) * prices[i + 1]
        total += 0.5 * (g0 + g1) * (k1 - k0)
    return total


def variance_swap_strike(S0, t, r,
                         put_strikes: Sequence[float], put_prices: Sequence[float],
                         call_strikes: Sequence[float], call_prices: Sequence[float],
                         split: float = None) -> float:
    """Fair variance-swap strike (annualized variance) by option replication.

    Args:
        S0: current spot.
        t: swap tenor in years.
        r: risk-free rate.
        put_strikes/put_prices: OTM puts, strikes strictly below ``split``.
        call_strikes/call_prices: OTM calls, strikes strictly above ``split``.
        split: the forward split level ``K*``. Defaults to the forward
            ``S0 e^{r t}``.

    Returns the fair strike as an annualized variance (multiply tenor and take
    sqrt for a vol number).
    """
    if t <= 0:
        raise ValueError("t must be positive")
    if split is None:
        split = S0 * math.exp(r * t)
    if split <= 0:
        raise ValueError("split level must be positive")

    growth = math.exp(r * t)

    # Strips weighted by 1/K^2 (the DDKZ replication weight).
    put_strip = _strip_integral(list(put_strikes), list(put_prices),
                                lambda K: 1.0 / (K * K)) if len(put_strikes) >= 2 else 0.0
    call_strip = _strip_integral(list(call_strikes), list(call_prices),
                                 lambda K: 1.0 / (K * K)) if len(call_strikes) >= 2 else 0.0

    # Log-contract decomposition around the split level K* (Demeterfi et al.).
    # K_var = (2/t)[ rt - (S0 e^{rt}/K* - 1) - ln(K*/S0) ] + (2 e^{rt}/t) * strips.
    forward_term = (S0 * growth / split) - 1.0 + math.log(split / S0)
    fair_var = (2.0 / t) * (r * t - forward_term) \
        + (2.0 * growth / t) * (put_strip + call_strip)
    return fair_var


def volatility_swap_strike(S0, t, r, put_strikes, put_prices,
                           call_strikes, call_prices, split=None) -> float:
    """Fair volatility-swap strike as ``sqrt(variance strike)``.

    This is the standard first-order proxy; it slightly overstates the true
    vol-swap strike because ``E[sqrt(var)] <= sqrt(E[var])`` (Jensen), the
    "convexity" or "vol-of-vol" adjustment, which requires a model to quantify.
    """
    var = variance_swap_strike(S0, t, r, put_strikes, put_prices,
                               call_strikes, call_prices, split)
    return math.sqrt(max(var, 0.0))


def volatility_swap_bounds_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401,
                                      width=8.0):
    """Bracket the fair volatility-swap strike from a smile: ``(lower, upper)``.

    The fair vol-swap strike ``E[sqrt(RV)]`` has no model-free replication (unlike
    the variance swap), but it is pinned between two computable levels:

      * upper -- ``sqrt(K_var)`` with ``K_var`` the variance-swap strike, since
        ``E[sqrt(RV)] <= sqrt(E[RV])`` by Jensen (the convexity/vol-of-vol gap).
      * lower -- the at-the-money-forward implied vol ``vol_fn(F)``. Carr-Lee
        (2009) show the fair vol-swap strike equals the ATMF implied vol exactly
        when spot/vol correlation is zero, and the ATMF vol lies below
        ``sqrt(K_var)`` whenever the smile is convex; it is the standard
        first-order proxy.

    Returns ``(atmf_vol, sqrt_var_strike)`` with ``lower <= upper``. A flat smile
    collapses the bracket to that flat vol (zero convexity). The true strike sits
    inside; the width ``upper - lower`` is the convexity premium the smile implies.
    """
    F = S0 * math.exp((r - q) * t)
    lower = vol_fn(F)
    var = variance_swap_from_smile(S0, t, r, vol_fn, q=q, n_strikes=n_strikes,
                                   width=width)
    upper = math.sqrt(max(var, 0.0))
    if lower > upper:
        # A non-convex (or numerically flat) smile can invert the pair; clamp so
        # the documented ordering always holds.
        lower = upper
    return lower, upper


def variance_swap_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401,
                             width=8.0, split=None):
    """Fair variance-swap strike from a volatility *smile* ``vol_fn(K)``.

    Builds the OTM option strip -- puts below the forward split, calls above --
    by pricing each strike at its smile vol ``vol_fn(K)`` with Black-Scholes,
    then feeds them to :func:`variance_swap_strike`. Convenient for marking a
    variance swap directly off a fitted smile (SVI, SABR, vanna-volga, ...).

    Args:
        vol_fn: callable ``vol_fn(K)`` returning the Black implied vol at strike.
        q: dividend yield (carry ``b = r - q``).
        n_strikes: number of strikes on each side; strikes span ``width``
            standard deviations of log-moneyness around the split.

    A flat smile returns exactly that flat variance (the model-free result).
    """
    from .bsm import call_price, put_price

    if t <= 0:
        raise ValueError("t must be positive")
    F = S0 * math.exp((r - q) * t)
    if split is None:
        split = F
    atm_vol = vol_fn(split)
    sd = atm_vol * math.sqrt(t)
    # Log-moneyness grid around the split, converted to strikes.
    put_ks, put_px, call_ks, call_px = [], [], [], []
    for i in range(n_strikes):
        x = -width * sd + 2.0 * width * sd * i / (n_strikes - 1)
        K = split * math.exp(x)
        v = vol_fn(K)
        if K < split:
            put_ks.append(K)
            put_px.append(put_price(S0, K, t, r, v, b=r - q))
        else:
            call_ks.append(K)
            call_px.append(call_price(S0, K, t, r, v, b=r - q))
    return variance_swap_strike(S0, t, r, put_ks, put_px, call_ks, call_px,
                                split=split)


def corridor_variance_swap_from_smile(S0, t, r, vol_fn, lower, upper, q=0.0,
                                      n_strikes=401, split=None):
    """Fair corridor variance-swap strike from a smile ``vol_fn(K)``.

    A corridor variance swap accrues realized variance only while the spot is in
    the corridor ``[lower, upper]``. By the Carr-Madan static-replication view
    this restricts the ``1/K^2``-weighted option strip to strikes inside the
    corridor (Carr & Lewis): the fair accrued variance is

        K_corr = (2 e^{r t} / t) * ( strip of OTM options with L <= K <= U ).

    ``vol_fn(K)`` prices each strip option with Black-Scholes at its smile vol.
    A corridor spanning the whole strip recovers (most of) the plain
    variance-swap strike; a narrower corridor accrues less variance.
    """
    from .bsm import call_price, put_price

    if t <= 0:
        raise ValueError("t must be positive")
    if not (0.0 < lower < upper):
        raise ValueError("need 0 < lower < upper")
    F = S0 * math.exp((r - q) * t)
    if split is None:
        split = F
    growth = math.exp(r * t)

    put_ks, put_px, call_ks, call_px = [], [], [], []
    # Uniform strike grid across the corridor.
    for i in range(n_strikes):
        K = lower + (upper - lower) * i / (n_strikes - 1)
        if K <= 0:
            continue
        v = vol_fn(K)
        if K < split:
            put_ks.append(K)
            put_px.append(put_price(S0, K, t, r, v, b=r - q))
        else:
            call_ks.append(K)
            call_px.append(call_price(S0, K, t, r, v, b=r - q))

    put_strip = _strip_integral(put_ks, put_px, lambda K: 1.0 / (K * K)) \
        if len(put_ks) >= 2 else 0.0
    call_strip = _strip_integral(call_ks, call_px, lambda K: 1.0 / (K * K)) \
        if len(call_ks) >= 2 else 0.0
    return (2.0 * growth / t) * (put_strip + call_strip)


def gamma_swap_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401, width=8.0,
                          split=None):
    """Fair gamma-swap strike from a smile ``vol_fn(K)``.

    A gamma (or "weighted variance") swap accrues ``(S_t / S0) d<ln S>`` -- each
    increment of realized variance weighted by the spot level -- so it is
    replicated by an option strip weighted ``1/K`` (the price-weighted version of
    the variance swap's ``1/K^2``), plus the matching log-contract terms
    (Carr-Lewis). Its fair strike is

        K_gamma = (2 e^{r t} / (S0 t)) * ( 1/K-weighted OTM strip )
                  + (2/t) * (r - q) * (e^{(r-q) t} - 1) / (r - q) ...  [drift term]

    Implemented from the price-weighted log contract; ``vol_fn(K)`` prices each
    option with Black-Scholes. A flat smile returns that flat variance.
    """
    from .bsm import call_price, put_price

    if t <= 0:
        raise ValueError("t must be positive")
    F = S0 * math.exp((r - q) * t)
    if split is None:
        split = F
    atm_vol = vol_fn(split)
    sd = atm_vol * math.sqrt(t)

    put_ks, put_px, call_ks, call_px = [], [], [], []
    for i in range(n_strikes):
        x = -width * sd + 2.0 * width * sd * i / (n_strikes - 1)
        K = split * math.exp(x)
        v = vol_fn(K)
        if K < split:
            put_ks.append(K)
            put_px.append(put_price(S0, K, t, r, v, b=r - q))
        else:
            call_ks.append(K)
            call_px.append(call_price(S0, K, t, r, v, b=r - q))

    # Price-weighted strip: weight 1/K (vs 1/K^2 for the plain variance swap).
    put_strip = _strip_integral(put_ks, put_px, lambda K: 1.0 / K) \
        if len(put_ks) >= 2 else 0.0
    call_strip = _strip_integral(call_ks, call_px, lambda K: 1.0 / K) \
        if len(call_ks) >= 2 else 0.0

    # Gamma-swap replication (Carr-Lewis): the price-weighted log contract has
    # no forward remainder (unlike the variance swap), so the fair accrued
    # variance is purely the 1/K-weighted strip scaled by 2 e^{rt}/(S0 t).
    growth = math.exp(r * t)
    return (2.0 * growth / (S0 * t)) * (put_strip + call_strip)


def forward_variance_swap_from_smile(S0, t1, t2, r, vol_fn1, vol_fn2, q=0.0,
                                     n_strikes=401, width=8.0):
    """Fair forward-start variance-swap strike over ``[t1, t2]`` from two smiles.

    Total (undiscounted) variance is additive in time, so the fair variance
    accrued between ``t1`` and ``t2`` is

        K_fwd = ( K_var(t2) * t2 - K_var(t1) * t1 ) / (t2 - t1),

    where ``K_var(t_i)`` is the spot-starting variance-swap strike replicated
    from the expiry-``t_i`` smile ``vol_fn_i(K)``. Requires ``0 <= t1 < t2``. A
    flat term structure of flat smiles returns that flat variance.
    """
    if not (0.0 <= t1 < t2):
        raise ValueError("need 0 <= t1 < t2")
    kv2 = variance_swap_from_smile(S0, t2, r, vol_fn2, q=q,
                                   n_strikes=n_strikes, width=width)
    if t1 == 0.0:
        return kv2
    kv1 = variance_swap_from_smile(S0, t1, r, vol_fn1, q=q,
                                   n_strikes=n_strikes, width=width)
    return (kv2 * t2 - kv1 * t1) / (t2 - t1)


def variance_term_structure(S0, r, expiries, vol_fns, q=0.0, n_strikes=401,
                            width=8.0):
    """Term structure of variance-swap strikes and the forward-variance curve.

    Args:
        expiries: increasing list of expiries.
        vol_fns: one smile ``vol_fn(K)`` per expiry (matching order).

    Returns ``(spot_var, forward_var)`` where ``spot_var[i]`` is the
    spot-starting variance-swap strike to ``expiries[i]`` and ``forward_var[i]``
    the annualized *forward* variance over ``(expiries[i-1], expiries[i]]``
    (``forward_var[0]`` = ``spot_var[0]``). By total-variance additivity
    ``forward_var[i] = (K_i t_i - K_{i-1} t_{i-1}) / (t_i - t_{i-1})``.

    A flat term structure of flat smiles gives a constant curve; a rising
    variance term structure gives positive, increasing forward variances.
    """
    if len(expiries) != len(vol_fns):
        raise ValueError("expiries and vol_fns must match in length")
    ts = list(expiries)
    if any(ts[i] >= ts[i + 1] for i in range(len(ts) - 1)):
        raise ValueError("expiries must be strictly increasing")

    spot_var = [variance_swap_from_smile(S0, ti, r, fn, q=q, n_strikes=n_strikes,
                                         width=width)
                for ti, fn in zip(ts, vol_fns)]
    forward_var = [spot_var[0]]
    for i in range(1, len(ts)):
        fv = (spot_var[i] * ts[i] - spot_var[i - 1] * ts[i - 1]) \
            / (ts[i] - ts[i - 1])
        forward_var.append(fv)
    return spot_var, forward_var
