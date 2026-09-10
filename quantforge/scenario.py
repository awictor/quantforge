"""Scenario / stress analysis: revalue a book across a grid of market shocks.

Traders and risk managers want to see the mark-to-market P&L of a position as
the underlying and its volatility move together — a "risk matrix" or stress
grid. This module reprices every contract in a :class:`Book` under each
(spot-shock, vol-shock) pair and returns the P&L surface, plus helpers to find
the worst-case scenario and the grid's delta/gamma profile.

All pricing goes through the exact scalar engine, so results agree with the
rest of the library. Shocks can be applied multiplicatively (relative) or as
absolute additive moves.
"""

import math
from dataclasses import dataclass, field
from typing import Sequence

from .portfolio import Contract, price_book


@dataclass(frozen=True)
class ScenarioGrid:
    spot_shocks: tuple          # the spot shocks applied (as given)
    vol_shocks: tuple           # the vol shocks applied (as given)
    pnl: tuple                  # pnl[i][j] for spot_shocks[i], vol_shocks[j]
    base_value: float           # book market value with no shock
    relative: bool

    def worst_case(self):
        """Return (spot_shock, vol_shock, pnl) of the largest loss on the grid."""
        worst = None
        for i, ss in enumerate(self.spot_shocks):
            for j, vs in enumerate(self.vol_shocks):
                p = self.pnl[i][j]
                if worst is None or p < worst[2]:
                    worst = (ss, vs, p)
        return worst

    def best_case(self):
        best = None
        for i, ss in enumerate(self.spot_shocks):
            for j, vs in enumerate(self.vol_shocks):
                p = self.pnl[i][j]
                if best is None or p > best[2]:
                    best = (ss, vs, p)
        return best

    def as_rows(self):
        """Yield (spot_shock, {vol_shock: pnl}) rows for tabular display."""
        for i, ss in enumerate(self.spot_shocks):
            yield ss, {vs: self.pnl[i][j] for j, vs in enumerate(self.vol_shocks)}


def _shock_contract(c: Contract, spot_shock, vol_shock, relative) -> Contract:
    if relative:
        new_S = c.S * (1.0 + spot_shock)
        new_sigma = c.sigma * (1.0 + vol_shock)
    else:
        new_S = c.S + spot_shock
        new_sigma = c.sigma + vol_shock
    new_S = max(new_S, 1e-12)
    new_sigma = max(new_sigma, 0.0)
    return Contract(S=new_S, K=c.K, t=c.t, r=c.r, sigma=new_sigma,
                    option_type=c.option_type, b=c.b, qty=c.qty,
                    multiplier=c.multiplier, label=c.label)


def stress_grid(contracts: Sequence[Contract], spot_shocks, vol_shocks,
                relative=True) -> ScenarioGrid:
    """Reprice a book across a Cartesian grid of spot and vol shocks.

    Args:
        contracts: the positions (signed qty, multiplier as in ``price_book``).
        spot_shocks: iterable of shocks to the underlying (e.g. [-0.1, 0, 0.1]).
        vol_shocks: iterable of shocks to volatility.
        relative: if True shocks are fractional (0.1 = +10%); if False they are
            absolute additive moves (spot in price units, vol in vol points).

    Returns a :class:`ScenarioGrid` whose ``pnl[i][j]`` is the change in book
    market value under ``(spot_shocks[i], vol_shocks[j])``.
    """
    contracts = [c.normalized() for c in contracts]
    spot_shocks = tuple(spot_shocks)
    vol_shocks = tuple(vol_shocks)
    base_value = price_book(contracts).net.market_value

    pnl = []
    for ss in spot_shocks:
        row = []
        for vs in vol_shocks:
            shocked = [_shock_contract(c, ss, vs, relative) for c in contracts]
            mv = price_book(shocked).net.market_value
            row.append(mv - base_value)
        pnl.append(tuple(row))
    return ScenarioGrid(spot_shocks=spot_shocks, vol_shocks=vol_shocks,
                        pnl=tuple(pnl), base_value=base_value, relative=relative)


def spot_ladder(contracts: Sequence[Contract], spot_shocks, relative=True):
    """A 1-D price ladder: book P&L vs spot shock only (vol unchanged).

    Returns a list of (spot_shock, pnl) pairs.
    """
    grid = stress_grid(contracts, spot_shocks, [0.0], relative=relative)
    return [(ss, grid.pnl[i][0]) for i, ss in enumerate(grid.spot_shocks)]
