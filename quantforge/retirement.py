"""Retirement drawdown: sustainable withdrawals, depletion, and glide paths.

Deterministic (constant-return) tools for the decumulation phase: how long a
portfolio lasts under a fixed or inflation-indexed withdrawal, the largest
sustainable withdrawal over a horizon, and a linear equity glide path. These are
the closed-form annuity relationships underlying rules of thumb like the "4% rule".
Pure standard library.
"""

import math


def portfolio_depletion_years(balance, annual_withdrawal, real_return):
    """Years a portfolio lasts under a constant real withdrawal.

    Solves the annuity exhaustion ``balance = W * (1 - (1+g)^{-n}) / g`` for ``n``
    with real growth ``g = real_return``:

        n = -ln(1 - g * balance / W) / ln(1 + g).

    Returns ``inf`` when the withdrawal is at or below the interest earned
    (``W <= g * balance``); at ``g = 0`` it is simply ``balance / W``.
    """
    if balance <= 0 or annual_withdrawal <= 0:
        raise ValueError("balance and annual_withdrawal must be positive")
    if real_return <= -1.0:
        raise ValueError("real_return must exceed -100%")
    g = real_return
    if abs(g) < 1e-12:
        return balance / annual_withdrawal
    if annual_withdrawal <= g * balance:
        return float("inf")
    return -math.log(1.0 - g * balance / annual_withdrawal) / math.log(1.0 + g)


def sustainable_withdrawal(balance, real_return, years):
    """Largest constant real withdrawal that exactly depletes over ``years``.

    The annuity payment ``W = balance * g / (1 - (1+g)^{-n})`` (``= balance / n`` at
    zero real return). The inverse of :func:`portfolio_depletion_years`.
    """
    if balance <= 0 or years <= 0:
        raise ValueError("balance and years must be positive")
    if real_return <= -1.0:
        raise ValueError("real_return must exceed -100%")
    g = real_return
    if abs(g) < 1e-12:
        return balance / years
    return balance * g / (1.0 - (1.0 + g) ** (-years))


def withdrawal_balance_path(balance, annual_withdrawal, nominal_return,
                            inflation, years):
    """Year-end balances under an inflation-indexed withdrawal.

    The withdrawal grows with ``inflation`` each year while the portfolio grows at
    the ``nominal_return``; withdrawals happen at year start. Returns the list of
    year-end balances (clipped at zero once depleted).
    """
    if balance < 0:
        raise ValueError("balance must be non-negative")
    if nominal_return <= -1.0 or inflation <= -1.0:
        raise ValueError("rates must exceed -100%")
    bal = balance
    w = annual_withdrawal
    out = []
    for _ in range(int(years)):
        bal = max(bal - w, 0.0) * (1.0 + nominal_return)
        out.append(bal)
        w *= 1.0 + inflation
    return out


def withdrawal_stream_pv(annual_withdrawal, real_discount_rate, years,
                         growth=0.0):
    """Present value of a (possibly growing) real withdrawal stream.

    Discounts ``years`` annual withdrawals -- ``annual_withdrawal`` growing at
    ``growth`` per year -- at ``real_discount_rate``, withdrawals at year start:

        PV = sum_{k=0}^{n-1} W (1+growth)^k / (1+r)^k.

    The capital needed to fund the stream. Rises with the withdrawal, the horizon,
    and the growth rate; falls with the discount rate.
    """
    if annual_withdrawal < 0:
        raise ValueError("annual_withdrawal must be non-negative")
    if real_discount_rate <= -1.0 or growth <= -1.0:
        raise ValueError("rates must exceed -100%")
    ratio = (1.0 + growth) / (1.0 + real_discount_rate)
    n = int(years)
    if abs(ratio - 1.0) < 1e-12:
        return annual_withdrawal * n
    return annual_withdrawal * (1.0 - ratio ** n) / (1.0 - ratio)


def ruin_probability_mc(balance, annual_withdrawal, mean_return, vol, years,
                        inflation=0.0, n_paths=10000, seed=20260911):
    """Probability of portfolio ruin under lognormal returns (Monte Carlo).

    Simulates annual real returns as ``N(mean_return, vol^2)``, taking an
    inflation-indexed withdrawal at the start of each year, and reports the
    fraction of paths that hit zero before ``years``. Rises with the withdrawal
    rate and volatility; near zero for withdrawals well below the mean return.
    Deterministic per seed.
    """
    if balance <= 0 or annual_withdrawal < 0:
        raise ValueError("balance must be positive, withdrawal non-negative")
    if vol < 0:
        raise ValueError("vol must be non-negative")
    if n_paths < 1 or years < 1:
        raise ValueError("n_paths and years must be positive")
    state = seed & 0xFFFFFFFF

    def _normal():
        nonlocal state
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        u1 = (state + 0.5) / 0x80000000
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        u2 = (state + 0.5) / 0x80000000
        return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)

    ruined = 0
    n = int(years)
    for _ in range(n_paths):
        bal = balance
        w = annual_withdrawal
        for _y in range(n):
            bal -= w
            if bal <= 0.0:
                ruined += 1
                break
            bal *= 1.0 + mean_return + vol * _normal()
            w *= 1.0 + inflation
    return ruined / n_paths


def glide_path_equity_weight(years_to_target, glide_years, start_equity,
                             end_equity):
    """Linear equity weight along a target-date glide path.

    Interpolates the equity allocation from ``start_equity`` (``glide_years`` out)
    down to ``end_equity`` at the target, clamped outside the glide window:

        w = end + (start - end) * clamp(years_to_target / glide_years, 0, 1).

    Declines monotonically as the target approaches.
    """
    if glide_years <= 0:
        raise ValueError("glide_years must be positive")
    if not (0.0 <= start_equity <= 1.0 and 0.0 <= end_equity <= 1.0):
        raise ValueError("equity weights must be in [0, 1]")
    frac = years_to_target / glide_years
    frac = min(max(frac, 0.0), 1.0)
    return end_equity + (start_equity - end_equity) * frac
