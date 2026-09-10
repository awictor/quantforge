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
