"""Merton (1976) jump-diffusion option pricing.

Merton adds lognormally-distributed jumps to the Black-Scholes diffusion: the
underlying jumps at Poisson rate ``lam`` per year, and each jump multiplies the
price by a lognormal factor with log-mean ``mu_j`` and log-std ``sigma_j``.

Conditional on exactly ``n`` jumps occurring before expiry, the terminal price
is again lognormal, so the option value is a Poisson-weighted sum of
Black-Scholes prices with adjusted volatility and carry:

    price = sum_{n=0..inf} e^{-lam' t} (lam' t)^n / n! * BSM(S, K, t, r_n, sigma_n)

where ``lam' = lam * (1 + k)``, ``k = exp(mu_j + 0.5 sigma_j^2) - 1`` is the
mean proportional jump size, and the n-jump diffusion parameters are

    sigma_n^2 = sigma^2 + n * sigma_j^2 / t
    b_n       = b - lam * k + n * (mu_j + 0.5 sigma_j^2) / t.

The series converges fast (terms decay like a Poisson tail); we sum until the
Poisson weight is negligible.
"""

import math

from .bsm import price as bsm_price, OptionType, _coerce_type, _validate


def merton_jump_price(S, K, t, r, sigma, lam, mu_j, sigma_j,
                      option_type=OptionType.CALL, b=None,
                      max_terms=200, tol=1e-12) -> float:
    """Price a European option under the Merton jump-diffusion model.

    Args:
        sigma: diffusion volatility (the continuous part).
        lam: jump intensity (expected number of jumps per year, >= 0).
        mu_j: mean of the log jump size.
        sigma_j: standard deviation of the log jump size (>= 0).
        b: cost of carry (defaults to r). The drift is compensated so the
            discounted asset is a martingale under the given carry.

    Returns the option price as the Poisson-weighted BSM series.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if lam < 0 or sigma_j < 0:
        raise ValueError("lam and sigma_j must be non-negative")
    if b is None:
        b = r

    if t == 0:
        return max(S - K, 0.0) if ot is OptionType.CALL else max(K - S, 0.0)

    # Mean proportional jump size. ln(1+k) = mu_j + 0.5 sigma_j^2.
    k = math.exp(mu_j + 0.5 * sigma_j * sigma_j) - 1.0
    ln1k = mu_j + 0.5 * sigma_j * sigma_j

    # Condition on n jumps (Poisson rate lam). The n-jump carry drifts by the
    # continuous compensation -lam*k plus n jumps of ln(1+k); Poisson-weighting
    # this by rate lam reproduces the martingale forward S*e^{b t} exactly.
    total = 0.0
    weight = math.exp(-lam * t)   # w_0 = e^{-lam t}
    for n in range(max_terms):
        var_n = sigma * sigma + n * sigma_j * sigma_j / t
        sigma_n = math.sqrt(var_n)
        b_n = b - lam * k + n * ln1k / t
        term = weight * bsm_price(S, K, t, r, sigma_n, ot, b=b_n)
        total += term
        # w_{n+1} = w_n * lam t / (n+1).
        weight *= (lam * t) / (n + 1)
        if n > lam * t and weight < tol:
            break

    return total


def merton_jump_greeks(S, K, t, r, sigma, lam, mu_j, sigma_j,
                       option_type=OptionType.CALL, b=None):
    """Greeks of a Merton jump-diffusion option by central finite differences.

    Differentiates :func:`merton_jump_price` for ``delta`` (dV/dS), ``gamma``
    (d2V/dS2), ``vega`` (dV/dsigma, the *diffusion*-vol sensitivity), and
    ``theta`` (calendar decay). At ``lam = 0`` (no jumps) the Greeks reduce to the
    vanilla Black-Scholes Greeks. Returns a dict with ``price`` and those fields.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if lam < 0 or sigma_j < 0:
        raise ValueError("lam and sigma_j must be non-negative")
    if b is None:
        b = r

    def px(S_=S, t_=t, sigma_=sigma):
        return merton_jump_price(S_, K, t_, r, sigma_, lam, mu_j, sigma_j, ot,
                                 b=b)

    base = px()
    hS = 1e-4 * S
    up, dn = px(S_=S + hS), px(S_=S - hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    hv = 1e-4
    vega = (px(sigma_=sigma + hv) - px(sigma_=sigma - hv)) / (2.0 * hv)
    ht = min(1e-4, 0.25 * t)
    theta = -(px(t_=t + ht) - px(t_=t - ht)) / (2.0 * ht)
    return {"price": base, "delta": delta, "gamma": gamma, "vega": vega,
            "theta": theta}


def merton_smile(S, strikes, t, r, sigma, lam, mu_j, sigma_j, b=None):
    """The Black-Scholes implied-vol smile a Merton jump-diffusion produces.

    Prices a European call at each strike under the jump-diffusion, then inverts
    each price to its Black-Scholes implied volatility, returning
    ``(log_moneyness, vol)`` pairs sorted by strike (log-moneyness on the forward
    ``F = S e^{b t}``). Jumps fatten the tails, so the smile curves up in the
    wings; a negative mean jump ``mu_j`` tilts it into a downward skew.
    """
    from .implied import implied_volatility
    if b is None:
        b = r
    F = S * math.exp(b * t)
    out = []
    for K in sorted(strikes):
        c = merton_jump_price(S, K, t, r, sigma, lam, mu_j, sigma_j,
                              OptionType.CALL, b=b)
        try:
            iv = implied_volatility(c, S, K, t, r, OptionType.CALL, b=b)
        except ValueError:
            continue
        out.append((math.log(K / F), iv))
    return out
