"""Vega term-structure bucketing for a multi-expiry option book.

Total book vega hides where the vol risk sits along the curve. This groups each
position's position-scaled vega into maturity buckets, so a desk sees whether it
is long front-month vol and short back-month, etc. Buckets are defined by their
upper-edge tenors; each contract lands in the first bucket whose edge is >= its
expiry (anything past the last edge goes to a final "beyond" bucket).

Works off :class:`~quantforge.portfolio.Contract` legs and the analytic vega, so
it composes with the rest of the engine.
"""

from bisect import bisect_left
from dataclasses import dataclass
from typing import Dict, List, Sequence

from .portfolio import Contract
from .bsm import vega as bs_vega


@dataclass(frozen=True)
class VegaBuckets:
    buckets: Dict[str, float]     # label -> net position-scaled vega
    total: float                  # sum over buckets (= net book vega)

    def as_list(self) -> List[tuple]:
        return list(self.buckets.items())


def _bucket_label(edges, i):
    if i == 0:
        return f"<={edges[0]:g}y"
    if i < len(edges):
        return f"{edges[i-1]:g}-{edges[i]:g}y"
    return f">{edges[-1]:g}y"


def vega_buckets(contracts: Sequence[Contract],
                 edges: Sequence[float] = (0.25, 0.5, 1.0, 2.0, 5.0)) -> VegaBuckets:
    """Bucket a book's position-scaled vega by expiry.

    Args:
        contracts: the book's legs (signed qty, multiplier as in ``price_book``).
        edges: sorted upper-edge tenors in years. A contract with expiry ``t``
            falls in the first bucket whose edge is >= ``t``; longer expiries go
            to the final ">last" bucket.

    Returns a :class:`VegaBuckets`. Vega is per 1.0 change in vol (divide by 100
    for per-vol-point), scaled by ``qty * multiplier``.
    """
    edges = list(edges)
    for i in range(len(edges) - 1):
        if edges[i + 1] <= edges[i]:
            raise ValueError("edges must be strictly increasing")

    labels = [_bucket_label(edges, i) for i in range(len(edges) + 1)]
    totals = {lab: 0.0 for lab in labels}

    for raw in contracts:
        c = raw.normalized()
        v = c.qty * c.multiplier * bs_vega(c.S, c.K, c.t, c.r, c.sigma, b=c.b)
        idx = bisect_left(edges, c.t)   # first edge >= t
        totals[labels[idx]] += v

    return VegaBuckets(buckets=totals, total=sum(totals.values()))
