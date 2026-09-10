"""Batch pricing and portfolio-level risk aggregation.

``price_book`` values a list of contracts in one call and returns per-contract
Greeks plus net book exposures. Still pure standard library — no NumPy.
"""

from dataclasses import dataclass, field
from typing import Iterable

from .bsm import greeks, Greeks, OptionType, _coerce_type


@dataclass(frozen=True)
class Contract:
    """A single option position.

    ``qty`` is signed: positive = long, negative = short. ``multiplier`` scales
    each contract to its notional (e.g. 100 for US equity options).
    """
    S: float
    K: float
    t: float
    r: float
    sigma: float
    option_type: OptionType = OptionType.CALL
    b: float = None
    qty: float = 1.0
    multiplier: float = 1.0
    label: str = ""

    def normalized(self) -> "Contract":
        return Contract(
            S=self.S, K=self.K, t=self.t, r=self.r, sigma=self.sigma,
            option_type=_coerce_type(self.option_type),
            b=self.b, qty=self.qty, multiplier=self.multiplier, label=self.label,
        )


@dataclass(frozen=True)
class Position:
    """A contract paired with its computed Greeks and position-scaled values."""
    contract: Contract
    greeks: Greeks

    @property
    def scale(self) -> float:
        return self.contract.qty * self.contract.multiplier

    @property
    def market_value(self) -> float:
        return self.scale * self.greeks.price

    @property
    def position_delta(self) -> float:
        return self.scale * self.greeks.delta

    @property
    def position_gamma(self) -> float:
        return self.scale * self.greeks.gamma

    @property
    def position_vega(self) -> float:
        return self.scale * self.greeks.vega

    @property
    def position_theta(self) -> float:
        return self.scale * self.greeks.theta

    @property
    def position_rho(self) -> float:
        return self.scale * self.greeks.rho


@dataclass(frozen=True)
class BookRisk:
    """Aggregate book-level exposures (position-scaled sums)."""
    market_value: float = 0.0
    delta: float = 0.0
    gamma: float = 0.0
    vega: float = 0.0
    theta: float = 0.0
    rho: float = 0.0


@dataclass
class Book:
    positions: list = field(default_factory=list)
    net: BookRisk = field(default_factory=BookRisk)


def price_book(contracts: Iterable[Contract]) -> Book:
    """Value every contract and aggregate net book Greeks.

    Returns a ``Book`` with per-position detail and a ``net`` ``BookRisk`` of
    position-scaled (qty * multiplier) sums.
    """
    positions = []
    mv = d = g = v = th = rh = 0.0
    for raw in contracts:
        c = raw.normalized()
        gk = greeks(c.S, c.K, c.t, c.r, c.sigma, c.option_type, c.b)
        pos = Position(contract=c, greeks=gk)
        positions.append(pos)
        mv += pos.market_value
        d += pos.position_delta
        g += pos.position_gamma
        v += pos.position_vega
        th += pos.position_theta
        rh += pos.position_rho
    return Book(
        positions=positions,
        net=BookRisk(market_value=mv, delta=d, gamma=g, vega=v, theta=th, rho=rh),
    )
