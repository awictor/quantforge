"""Commodity forward pricing under the cost-of-carry model.

A storable commodity's forward is the spot compounded at the financing rate plus
storage cost and net of the convenience yield -- the benefit of holding the
physical good rather than a claim on it:

    F(T) = S * exp((r + u - y) * T)

with ``r`` the risk-free rate, ``u`` the (continuous, proportional) storage cost,
and ``y`` the convenience yield. When ``y > r + u`` the curve is in backwardation
(forwards below spot); otherwise it is in contango. This module prices the
forward, inverts the market forward for the implied convenience yield, builds a
forward curve, and classifies the curve shape. Pure standard library.
"""

import math


def commodity_forward(spot, r, maturity, storage_cost=0.0, convenience_yield=0.0):
    """Cost-of-carry forward ``S * exp((r + u - y) * T)``.

    ``storage_cost`` (u) and ``convenience_yield`` (y) are continuous proportional
    rates. Storage lifts the forward (a cost of carrying the physical), the
    convenience yield lowers it (a benefit of holding it).
    """
    if spot <= 0:
        raise ValueError("spot must be positive")
    if maturity < 0:
        raise ValueError("maturity must be non-negative")
    return spot * math.exp((r + storage_cost - convenience_yield) * maturity)


def implied_convenience_yield(spot, forward, r, maturity, storage_cost=0.0):
    """Convenience yield implied by a market forward (inverts the carry formula).

    Solves ``F = S exp((r + u - y) T)`` for ``y``:
    ``y = r + u - ln(F/S)/T``. Inverse of :func:`commodity_forward`.
    """
    if spot <= 0 or forward <= 0:
        raise ValueError("spot and forward must be positive")
    if maturity <= 0:
        raise ValueError("maturity must be positive")
    return r + storage_cost - math.log(forward / spot) / maturity


def implied_storage_cost(spot, forward, r, maturity, convenience_yield=0.0):
    """Storage cost implied by a market forward (inverts the carry formula).

    ``u = ln(F/S)/T - r + y``. Inverse of :func:`commodity_forward` in ``u``.
    """
    if spot <= 0 or forward <= 0:
        raise ValueError("spot and forward must be positive")
    if maturity <= 0:
        raise ValueError("maturity must be positive")
    return math.log(forward / spot) / maturity - r + convenience_yield


def net_cost_of_carry(r, storage_cost=0.0, convenience_yield=0.0):
    """Net proportional carry rate ``r + u - y`` (the forward's growth rate)."""
    return r + storage_cost - convenience_yield


def commodity_forward_curve(spot, r, maturities, storage_cost=0.0,
                            convenience_yield=0.0):
    """Forward prices across a list of maturities under one carry rate.

    Returns ``[(T, F(T)), ...]`` from :func:`commodity_forward` at each maturity.
    """
    return [(T, commodity_forward(spot, r, T, storage_cost, convenience_yield))
            for T in maturities]


def schwartz_log_mean(spot, kappa, alpha_star, maturity):
    """Risk-neutral mean of log-spot under the Schwartz one-factor model.

    Log-spot ``X = ln S`` follows a mean-reverting Ornstein-Uhlenbeck process
    ``dX = kappa (alpha_star - X) dt + sigma dW`` under the pricing measure, where
    ``alpha_star`` is the risk-adjusted long-run log level. The conditional mean is

        E[X_T] = e^{-kappa T} ln S + (1 - e^{-kappa T}) alpha_star,

    starting at ``ln S`` for ``T = 0`` and relaxing to ``alpha_star`` as
    ``T -> inf``.
    """
    if spot <= 0:
        raise ValueError("spot must be positive")
    if kappa <= 0:
        raise ValueError("kappa must be positive")
    if maturity < 0:
        raise ValueError("maturity must be non-negative")
    decay = math.exp(-kappa * maturity)
    return decay * math.log(spot) + (1.0 - decay) * alpha_star


def schwartz_log_variance(sigma, kappa, maturity):
    """Variance of log-spot under the Schwartz one-factor model.

    ``Var[X_T] = sigma^2 (1 - e^{-2 kappa T}) / (2 kappa)`` -- zero at ``T = 0``,
    rising monotonically to the stationary ``sigma^2 / (2 kappa)`` as
    ``T -> inf``.
    """
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if kappa <= 0:
        raise ValueError("kappa must be positive")
    if maturity < 0:
        raise ValueError("maturity must be non-negative")
    return sigma * sigma * (1.0 - math.exp(-2.0 * kappa * maturity)) / (2.0 * kappa)


def schwartz_forward(spot, kappa, alpha_star, sigma, maturity):
    """Commodity forward under the Schwartz (1997) one-factor model.

    With log-spot lognormal, ``F(T) = exp(E[X_T] + 0.5 Var[X_T])`` from
    :func:`schwartz_log_mean` and :func:`schwartz_log_variance`. Equals the spot
    at ``T = 0`` and converges to the risk-neutral long-run forward
    ``exp(alpha_star + sigma^2/(4 kappa))`` as ``T -> inf`` -- the mean-reverting
    alternative to the constant-carry :func:`commodity_forward`.
    """
    mean = schwartz_log_mean(spot, kappa, alpha_star, maturity)
    var = schwartz_log_variance(sigma, kappa, maturity)
    return math.exp(mean + 0.5 * var)


def schwartz_option(spot, kappa, alpha_star, sigma, strike, r, expiry,
                    is_call=True):
    """European spot option under the Schwartz one-factor model.

    At expiry the spot is lognormal with mean :func:`schwartz_log_mean` and
    variance :func:`schwartz_log_variance`, so the option is a Black-style price
    off the model forward ``F* = schwartz_forward`` and total variance
    ``v = Var[X_T]``:

        d1 = (ln(F*/K) + 0.5 v) / sqrt(v),  d2 = d1 - sqrt(v)
        call = e^{-r T} [F* Phi(d1) - K Phi(d2)]
        put  = e^{-r T} [K Phi(-d2) - F* Phi(-d1)]

    Put and call satisfy ``C - P = e^{-r T} (F* - K)``. At zero variance the price
    is the discounted intrinsic on ``F*``.
    """
    from .mathfns import norm_cdf
    if strike <= 0:
        raise ValueError("strike must be positive")
    if expiry < 0:
        raise ValueError("expiry must be non-negative")
    fwd = schwartz_forward(spot, kappa, alpha_star, sigma, expiry)
    disc = math.exp(-r * expiry)
    var = schwartz_log_variance(sigma, kappa, expiry)
    if var <= 0.0:
        intrinsic = max(fwd - strike, 0.0) if is_call else max(strike - fwd, 0.0)
        return disc * intrinsic
    vsqrt = math.sqrt(var)
    d1 = (math.log(fwd / strike) + 0.5 * var) / vsqrt
    d2 = d1 - vsqrt
    if is_call:
        return disc * (fwd * norm_cdf(d1) - strike * norm_cdf(d2))
    return disc * (strike * norm_cdf(-d2) - fwd * norm_cdf(-d1))


def mean_reversion_half_life(kappa):
    """Half-life of mean reversion ``ln(2) / kappa`` (years).

    Time for a shock to log-spot to decay to half its size under the Schwartz OU
    dynamics; falls as the mean-reversion speed ``kappa`` rises.
    """
    if kappa <= 0:
        raise ValueError("kappa must be positive")
    return math.log(2.0) / kappa


def schwartz_implied_alpha(spot, forward, kappa, sigma, maturity):
    """Risk-neutral long-run log level implied by a single forward quote.

    Inverts :func:`schwartz_forward` for ``alpha_star``:

        alpha_star = [ln F - 0.5 Var[X_T] - e^{-kappa T} ln S] / (1 - e^{-kappa T}).

    Requires ``maturity > 0`` (at ``T = 0`` the forward carries no information
    about the long-run level). Round-trips with :func:`schwartz_forward`.
    """
    if spot <= 0 or forward <= 0:
        raise ValueError("spot and forward must be positive")
    if kappa <= 0:
        raise ValueError("kappa must be positive")
    if maturity <= 0:
        raise ValueError("maturity must be positive")
    decay = math.exp(-kappa * maturity)
    var = schwartz_log_variance(sigma, kappa, maturity)
    return (math.log(forward) - 0.5 * var - decay * math.log(spot)) / (1.0 - decay)


def margrabe_exchange_option(f1, f2, sigma1, sigma2, rho, r, expiry):
    """Margrabe (1978) option to exchange asset 2 for asset 1, on forwards.

    Exact closed form for the payoff ``max(F1 - F2, 0)`` (a zero-strike spread
    option). With the spread volatility
    ``sigma = sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)``:

        price = e^{-r T} [F1 Phi(d1) - F2 Phi(d2)]
        d1 = (ln(F1/F2) + 0.5 sigma^2 T) / (sigma sqrt(T)),  d2 = d1 - sigma sqrt(T)

    The building block the :func:`kirk_spread_option` approximation reduces to at
    zero strike.
    """
    from .mathfns import norm_cdf
    if f1 <= 0 or f2 <= 0:
        raise ValueError("forwards must be positive")
    if expiry < 0:
        raise ValueError("expiry must be non-negative")
    disc = math.exp(-r * expiry)
    var = sigma1 * sigma1 - 2.0 * rho * sigma1 * sigma2 + sigma2 * sigma2
    if var <= 0.0 or expiry == 0.0:
        return disc * max(f1 - f2, 0.0)
    vsqrt = math.sqrt(var * expiry)
    d1 = (math.log(f1 / f2) + 0.5 * var * expiry) / vsqrt
    d2 = d1 - vsqrt
    return disc * (f1 * norm_cdf(d1) - f2 * norm_cdf(d2))


def kirk_spread_option(f1, f2, strike, sigma1, sigma2, rho, r, expiry,
                       is_call=True):
    """Kirk (1995) approximation for a spread option on two forwards.

    Prices ``max(F1 - F2 - K, 0)`` (call) or ``max(K - (F1 - F2), 0)`` (put) by
    treating ``F2 + K`` as a single lognormal asset and applying Black with an
    effective spread volatility

        sigma_K = sqrt(sigma1^2 - 2 rho sigma1 sigma2 w + sigma2^2 w^2),
        w = F2 / (F2 + K).

    At ``K = 0`` (``w = 1``) it collapses exactly to the
    :func:`margrabe_exchange_option`. Call and put satisfy
    ``C - P = e^{-r T} (F1 - F2 - K)``.
    """
    from .mathfns import norm_cdf
    if f1 <= 0 or f2 <= 0:
        raise ValueError("forwards must be positive")
    if f2 + strike <= 0:
        raise ValueError("F2 + strike must be positive for the Kirk approximation")
    if expiry < 0:
        raise ValueError("expiry must be non-negative")
    disc = math.exp(-r * expiry)
    spread = f1 - f2 - strike
    denom = f2 + strike
    w = f2 / denom
    var = sigma1 * sigma1 - 2.0 * rho * sigma1 * sigma2 * w + sigma2 * sigma2 * w * w
    if var <= 0.0 or expiry == 0.0:
        intrinsic = max(spread, 0.0) if is_call else max(-spread, 0.0)
        return disc * intrinsic
    vsqrt = math.sqrt(var * expiry)
    d1 = (math.log(f1 / denom) + 0.5 * var * expiry) / vsqrt
    d2 = d1 - vsqrt
    if is_call:
        return disc * (f1 * norm_cdf(d1) - denom * norm_cdf(d2))
    return disc * (denom * norm_cdf(-d2) - f1 * norm_cdf(-d1))


def bachelier_spread_option(f1, f2, strike, sigma1, sigma2, rho, r, expiry,
                            is_call=True):
    """Bachelier (normal-model) spread option on two forwards.

    Models each forward as arithmetic Brownian motion, so the spread ``F1 - F2``
    is normal with volatility
    ``sigma = sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)`` (absolute, price
    units). Unlike the lognormal :func:`kirk_spread_option` this prices spreads
    that are or can go negative -- the norm for crack and location spreads. With
    ``m = F1 - F2 - K`` and ``s = sigma sqrt(T)``:

        call = e^{-r T} [m Phi(m/s) + s phi(m/s)]
        put  = e^{-r T} [-m Phi(-m/s) + s phi(m/s)]

    Put and call satisfy ``C - P = e^{-r T} (F1 - F2 - K)``.
    """
    from .mathfns import norm_cdf
    if expiry < 0:
        raise ValueError("expiry must be non-negative")
    if sigma1 < 0 or sigma2 < 0:
        raise ValueError("volatilities must be non-negative")
    disc = math.exp(-r * expiry)
    m = f1 - f2 - strike
    var = sigma1 * sigma1 - 2.0 * rho * sigma1 * sigma2 + sigma2 * sigma2
    if var <= 0.0 or expiry == 0.0:
        intrinsic = max(m, 0.0) if is_call else max(-m, 0.0)
        return disc * intrinsic
    s = math.sqrt(var * expiry)
    d = m / s
    pdf = math.exp(-0.5 * d * d) / math.sqrt(2.0 * math.pi)
    if is_call:
        return disc * (m * norm_cdf(d) + s * pdf)
    return disc * (-m * norm_cdf(-d) + s * pdf)


def spread_option_mc(f1, f2, strike, sigma1, sigma2, rho, r, expiry,
                     n_paths=100000, seed=12345, is_call=True):
    """Monte Carlo price of a spread option under bivariate lognormal forwards.

    Simulates ``F1 e^{-0.5 sigma1^2 T + sigma1 sqrt(T) Z1}`` and the analogous
    ``F2`` with correlated normals ``corr(Z1, Z2) = rho`` (Cholesky), averaging the
    discounted payoff ``max(F1 - F2 - K, 0)`` (call) or its put. An independent
    reference for the :func:`kirk_spread_option` approximation. Uses a
    deterministic linear-congruential stream so results are reproducible.
    """
    if f1 <= 0 or f2 <= 0:
        raise ValueError("forwards must be positive")
    if expiry < 0:
        raise ValueError("expiry must be non-negative")
    if n_paths < 1:
        raise ValueError("n_paths must be positive")
    # Deterministic LCG -> uniforms -> Box-Muller normals (no external deps).
    state = seed & 0xFFFFFFFF
    def _unif():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        return (state + 0.5) / 0x80000000
    disc = math.exp(-r * expiry)
    sq = math.sqrt(expiry)
    chol = math.sqrt(max(1.0 - rho * rho, 0.0))
    drift1 = -0.5 * sigma1 * sigma1 * expiry
    drift2 = -0.5 * sigma2 * sigma2 * expiry
    total = 0.0
    for _ in range(n_paths):
        u1 = _unif()
        u2 = _unif()
        rmag = math.sqrt(-2.0 * math.log(u1))
        z1 = rmag * math.cos(2.0 * math.pi * u2)
        z2 = rmag * math.sin(2.0 * math.pi * u2)
        w2 = rho * z1 + chol * z2
        s1 = f1 * math.exp(drift1 + sigma1 * sq * z1)
        s2 = f2 * math.exp(drift2 + sigma2 * sq * w2)
        spread = s1 - s2 - strike
        payoff = max(spread, 0.0) if is_call else max(-spread, 0.0)
        total += payoff
    return disc * total / n_paths


def roll_yield(near_forward, far_forward, t_near, t_far):
    """Annualized roll yield between two forwards ``ln(F_near/F_far)/(t_far-t_near)``.

    The return earned rolling a long position from the far to the near contract as
    time passes, assuming spot is unchanged. Positive in backwardation (near above
    far) and negative in contango, so its sign is the opposite of the
    :func:`commodity_calendar_spread` sign.
    """
    if near_forward <= 0 or far_forward <= 0:
        raise ValueError("forwards must be positive")
    if t_far <= t_near:
        raise ValueError("t_far must exceed t_near")
    return math.log(near_forward / far_forward) / (t_far - t_near)


def carry_roll_yield(r, storage_cost=0.0, convenience_yield=0.0):
    """Roll yield implied by the cost-of-carry model ``y - r - u`` (= -net carry).

    Under constant carry ``F(T) = S e^{(r+u-y)T}`` the roll yield is exactly the
    negative net carry :func:`net_cost_of_carry`, so it equals the convenience
    yield net of financing and storage. Positive precisely in backwardation.
    """
    return -net_cost_of_carry(r, storage_cost, convenience_yield)


def schwartz_futures_volatility(sigma, kappa, maturity):
    """Instantaneous return volatility of the ``maturity``-future under Schwartz.

    Because log-spot mean-reverts, the futures return volatility decays with time
    to maturity: ``sigma_F(T) = sigma e^{-kappa T}``. It equals the spot vol
    ``sigma`` for the front (``T = 0``) and falls for longer maturities -- the
    Samuelson effect (near contracts more volatile than deferred).
    """
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if kappa <= 0:
        raise ValueError("kappa must be positive")
    if maturity < 0:
        raise ValueError("maturity must be non-negative")
    return sigma * math.exp(-kappa * maturity)


def commodity_calendar_spread(spot, r, t_near, t_far, storage_cost=0.0,
                              convenience_yield=0.0):
    """Far-minus-near forward spread under one carry rate.

    ``F(t_far) - F(t_near)`` from :func:`commodity_forward`. Positive in contango
    (net carry ``r + u - y > 0``, the far contract richer) and negative in
    backwardation, so its sign matches :func:`net_cost_of_carry`.
    """
    if t_far <= t_near:
        raise ValueError("t_far must exceed t_near")
    near = commodity_forward(spot, r, t_near, storage_cost, convenience_yield)
    far = commodity_forward(spot, r, t_far, storage_cost, convenience_yield)
    return far - near


def convenience_yield_curve(spot, r, forward_quotes, storage_cost=0.0):
    """Per-tenor convenience yields implied by a forward strip.

    ``forward_quotes`` is ``[(T, F(T)), ...]``. Inverts each quote with
    :func:`implied_convenience_yield` at a common ``storage_cost``, so recomputing
    the forward at each ``(T, y_T)`` reprices the input strip exactly.
    """
    if spot <= 0:
        raise ValueError("spot must be positive")
    return [(T, implied_convenience_yield(spot, F, r, T, storage_cost))
            for T, F in forward_quotes]


def seasonal_forward(spot, r, maturity, seasonal_factor, storage_cost=0.0,
                     convenience_yield=0.0):
    """Cost-of-carry forward scaled by a multiplicative seasonal factor.

    ``seasonal_factor * commodity_forward(...)`` -- lifts or discounts the carry
    forward for the delivery month's seasonal pattern (e.g. gas in winter). With
    a :func:`quantforge.normalize_seasonal_factors` factor the annual average is
    unchanged.
    """
    if seasonal_factor <= 0:
        raise ValueError("seasonal_factor must be positive")
    return seasonal_factor * commodity_forward(spot, r, maturity, storage_cost,
                                               convenience_yield)


def is_backwardation(r, storage_cost=0.0, convenience_yield=0.0):
    """True when the curve is in backwardation (``y > r + u``, forwards below spot).

    Backwardation occurs when the convenience yield exceeds the financing-plus-
    storage carry, so the net carry :func:`net_cost_of_carry` is negative and
    forwards fall with maturity.
    """
    return net_cost_of_carry(r, storage_cost, convenience_yield) < 0.0
