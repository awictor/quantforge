"""Book-level aggregate Greeks, including second order.

:func:`quantforge.price_book` already gives net first-order Greeks. This adds
the position-scaled net *second-order* exposures a vol desk manages —
net vanna, vomma/volga, charm, veta, speed, zomma, and color — summed across a
book of :class:`~quantforge.portfolio.Contract` legs.

Each is scaled by ``qty * multiplier`` per leg, so the numbers are the book's
actual sensitivity, not per-unit.
"""

from dataclasses import dataclass
from typing import Iterable

from .portfolio import Contract
from .greeks2 import vanna, vomma, charm, veta, speed, zomma, color


@dataclass(frozen=True)
class BookSecondOrder:
    vanna: float = 0.0
    vomma: float = 0.0
    charm: float = 0.0
    veta: float = 0.0
    speed: float = 0.0
    zomma: float = 0.0
    color: float = 0.0


def book_second_order(contracts: Iterable[Contract]) -> BookSecondOrder:
    """Aggregate position-scaled second-order Greeks across a book.

    charm/veta/color are calendar-convention (per year), matching the scalar
    functions in :mod:`quantforge.greeks2`.
    """
    va = vo = ch = ve = sp = zo = co = 0.0
    for raw in contracts:
        c = raw.normalized()
        scale = c.qty * c.multiplier
        va += scale * vanna(c.S, c.K, c.t, c.r, c.sigma, b=c.b)
        vo += scale * vomma(c.S, c.K, c.t, c.r, c.sigma, b=c.b)
        ch += scale * charm(c.S, c.K, c.t, c.r, c.sigma, c.option_type, b=c.b)
        ve += scale * veta(c.S, c.K, c.t, c.r, c.sigma, b=c.b)
        sp += scale * speed(c.S, c.K, c.t, c.r, c.sigma, b=c.b)
        zo += scale * zomma(c.S, c.K, c.t, c.r, c.sigma, b=c.b)
        co += scale * color(c.S, c.K, c.t, c.r, c.sigma, b=c.b)
    return BookSecondOrder(vanna=va, vomma=vo, charm=ch, veta=ve,
                           speed=sp, zomma=zo, color=co)
