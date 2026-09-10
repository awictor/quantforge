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


@dataclass(frozen=True)
class ThetaCarry:
    theta: float            # net position-scaled theta (per year, calendar)
    gamma_rent: float       # -0.5 * Gamma * sigma^2 * S^2  (theta from convexity)
    residual: float         # theta - gamma_rent (drift / financing carry)


def theta_carry_report(contracts):
    """Decompose a book's net theta into gamma-rent and a residual carry term.

    The theta-gamma relationship says an option's time decay is dominated by the
    "gamma rent" paid on convexity: for a single underlying with volatility
    ``sigma``, ``theta ~= -0.5 * Gamma * sigma^2 * S^2`` plus a smaller
    drift/financing residual (rho- and dividend-carry effects).

    Sums the position-scaled net theta and the gamma-rent term across the book;
    the residual is the difference. When every leg shares one spot/vol (a
    single-name book) the gamma-rent uses that common S and sigma.
    """
    from .bsm import theta as bs_theta, gamma as bs_gamma

    net_theta = 0.0
    gamma_rent = 0.0
    for raw in contracts:
        c = raw.normalized()
        scale = c.qty * c.multiplier
        net_theta += scale * bs_theta(c.S, c.K, c.t, c.r, c.sigma, c.option_type, b=c.b)
        g = bs_gamma(c.S, c.K, c.t, c.r, c.sigma, b=c.b)
        gamma_rent += scale * (-0.5 * g * c.sigma * c.sigma * c.S * c.S)
    return ThetaCarry(theta=net_theta, gamma_rent=gamma_rent,
                      residual=net_theta - gamma_rent)


@dataclass(frozen=True)
class BumpGreeks:
    price: float
    delta: float
    gamma: float
    vega: float
    theta: float            # calendar (per year)


def book_bump_greeks(contracts, dS_frac=1e-3, dvol=1e-4, dt=1e-4):
    """Net book Greeks by bumping the shared market and repricing (model-free).

    Applies a common shock to every leg's spot, volatility, and time-to-expiry
    and reprices the whole book via :func:`quantforge.price_book`, so the net
    delta/gamma/vega/theta come out numerically without needing analytic Greeks
    for each instrument. Assumes all legs share one underlying and vol (a
    single-name book), the usual case for this kind of check.

    ``dS_frac`` is the relative spot bump; ``dvol`` and ``dt`` are absolute.
    Returns a :class:`BumpGreeks`.
    """
    from .portfolio import Contract, price_book

    contracts = [c.normalized() for c in contracts]
    if not contracts:
        raise ValueError("book is empty")
    S0 = contracts[0].S

    def reprice(dS=0.0, dv=0.0, dtau=0.0):
        legs = [Contract(S=c.S + dS, K=c.K, t=max(c.t - dtau, 1e-9), r=c.r,
                         sigma=max(c.sigma + dv, 1e-9), option_type=c.option_type,
                         b=c.b, qty=c.qty, multiplier=c.multiplier, label=c.label)
                for c in contracts]
        return price_book(legs).net.market_value

    base = reprice()
    hS = dS_frac * S0
    up, dn = reprice(dS=hS), reprice(dS=-hS)
    delta = (up - dn) / (2.0 * hS)
    gamma = (up - 2.0 * base + dn) / (hS * hS)
    vega = (reprice(dv=dvol) - reprice(dv=-dvol)) / (2.0 * dvol)
    # Calendar theta: value change as time advances (t decreases).
    theta = (reprice(dtau=dt) - base) / dt
    return BumpGreeks(price=base, delta=delta, gamma=gamma, vega=vega, theta=theta)
