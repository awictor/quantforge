"""Position-sizing helpers for hedging a book to a target Greek.

Given a :class:`~quantforge.portfolio.Book` and a hedging instrument (an
option or the underlying), these solve the quantity of the instrument needed to
drive a chosen net Greek to a target (usually zero):

  * ``delta_hedge_shares`` — shares of the underlying to zero net delta.
  * ``vega_neutral_quantity`` — units of a hedge option to zero net vega.
  * ``gamma_neutral_quantity`` — units of a hedge option to zero net gamma.
  * ``neutralize`` — quantity of a hedge option to hit a target for any of
    delta / gamma / vega.

All work off the position-scaled net Greeks already computed by
``price_book`` and the per-unit Greeks of the hedge, so they compose with the
rest of the engine.
"""

from .bsm import greeks, OptionType, _coerce_type
from .portfolio import Book


def _hedge_greeks(S, K, t, r, sigma, option_type, b):
    g = greeks(S, K, t, r, sigma, option_type, b)
    return g


def delta_hedge_shares(book: Book) -> float:
    """Shares of the underlying to add to zero the book's net delta.

    A share has delta 1, so the hedge is ``-net_delta`` shares (negative = sell).
    """
    return -book.net.delta


def neutralize(book: Book, greek: str, S, K, t, r, sigma,
               option_type=OptionType.CALL, b=None, multiplier=1.0,
               target: float = 0.0) -> float:
    """Units of the hedge option to move ``greek`` to ``target`` (default 0).

    ``greek`` is one of "delta", "gamma", "vega". Returns the signed quantity
    (in option units, before the multiplier is applied to notionals): solving
    ``net_greek + qty * multiplier * hedge_greek = target``.
    """
    if greek not in ("delta", "gamma", "vega"):
        raise ValueError("greek must be one of delta, gamma, vega")
    ot = _coerce_type(option_type)
    if b is None:
        b = r
    g = _hedge_greeks(S, K, t, r, sigma, ot, b)
    per_unit = {"delta": g.delta, "gamma": g.gamma, "vega": g.vega}[greek]
    if abs(per_unit) < 1e-15:
        raise ValueError(f"hedge instrument has ~zero {greek}; cannot neutralize")
    net = getattr(book.net, greek)
    return (target - net) / (multiplier * per_unit)


def vega_neutral_quantity(book: Book, S, K, t, r, sigma,
                          option_type=OptionType.CALL, b=None,
                          multiplier=1.0) -> float:
    """Units of a hedge option that zero the book's net vega."""
    return neutralize(book, "vega", S, K, t, r, sigma, option_type, b, multiplier)


def gamma_neutral_quantity(book: Book, S, K, t, r, sigma,
                           option_type=OptionType.CALL, b=None,
                           multiplier=1.0) -> float:
    """Units of a hedge option that zero the book's net gamma."""
    return neutralize(book, "gamma", S, K, t, r, sigma, option_type, b, multiplier)


def kelly_fraction_binary(win_prob, win_payoff, loss_amount=1.0):
    """Kelly fraction for a binary bet.

    Args:
        win_prob: probability of winning, p in (0, 1).
        win_payoff: net amount won per unit staked on a win (the "b" in b-to-1
            odds).
        loss_amount: amount lost per unit staked on a loss (default 1).

    Returns the fraction of bankroll to wager: ``f = (p*b - q*loss) / (b*loss)``
    where ``q = 1 - p``. A non-positive result means the bet has no edge; the
    optimal stake is then zero (returned as a negative/zero fraction for the
    caller to clamp).
    """
    if not (0.0 < win_prob < 1.0):
        raise ValueError("win_prob must be in (0, 1)")
    if win_payoff <= 0 or loss_amount <= 0:
        raise ValueError("payoffs must be positive")
    q = 1.0 - win_prob
    return (win_prob * win_payoff - q * loss_amount) / (win_payoff * loss_amount)


def kelly_fraction_continuous(expected_excess_return, variance, fraction=1.0):
    """Continuous Kelly allocation for a normally-distributed return.

    For a return with mean excess ``mu`` (over the risk-free rate) and variance
    ``sigma^2``, the growth-optimal leverage is ``f* = mu / sigma^2``. Multiply
    by ``fraction`` for fractional Kelly (e.g. 0.5 for half-Kelly, which trades
    a little growth for much lower drawdown).
    """
    if variance <= 0:
        raise ValueError("variance must be positive")
    return fraction * expected_excess_return / variance


def kelly_growth_rate(expected_excess_return, variance, leverage):
    """Expected log-growth rate at a given leverage (continuous Kelly).

    ``g(f) = f*mu - 0.5 * f^2 * sigma^2``. Maximized at the full-Kelly leverage
    ``f* = mu / sigma^2``; used to compare fractional-Kelly choices.
    """
    return leverage * expected_excess_return - 0.5 * leverage * leverage * variance


def kelly_fractions_multivariate(mean_excess_returns, cov, fraction=1.0):
    """Growth-optimal Kelly allocation across correlated assets.

    For a vector of excess returns with mean ``mu`` and covariance ``Sigma``, the
    continuous multivariate Kelly criterion maximizes the expected log-growth
    ``f . mu - 0.5 f . Sigma f``; the optimum is ``f* = Sigma^{-1} mu``. Returns
    the leverage vector, scaled by ``fraction`` for fractional Kelly.

    Reduces to the scalar ``mu / sigma^2`` for a single asset, and to the
    per-asset Kelly fractions when the covariance is diagonal (uncorrelated
    assets). Requires a positive-definite ``Sigma``.
    """
    from .portopt import _invert, _matvec

    n = len(mean_excess_returns)
    if len(cov) != n or any(len(row) != n for row in cov):
        raise ValueError("cov must be square and match mean_excess_returns")
    inv = _invert(cov)
    f = _matvec(inv, list(map(float, mean_excess_returns)))
    return [fraction * fi for fi in f]


def kelly_growth_rate_multivariate(mean_excess_returns, cov, leverages):
    """Expected log-growth rate of a multivariate allocation.

    ``g(f) = f . mu - 0.5 f . Sigma f``. Maximized at
    ``f* = Sigma^{-1} mu`` (:func:`kelly_fractions_multivariate` with
    ``fraction=1``); used to compare fractional-Kelly leverage vectors.
    """
    from .portopt import _matvec

    mu = list(map(float, mean_excess_returns))
    f = list(map(float, leverages))
    if len(f) != len(mu):
        raise ValueError("leverages must match mean_excess_returns")
    linear = sum(f[i] * mu[i] for i in range(len(f)))
    Sf = _matvec(cov, f)
    quad = sum(f[i] * Sf[i] for i in range(len(f)))
    return linear - 0.5 * quad
