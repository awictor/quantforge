"""Option-strategy builders: standard multi-leg structures with net risk.

Each builder returns a :class:`~quantforge.portfolio.Book` of the constituent
legs, so the existing engine gives net price and Greeks for free
(:func:`quantforge.price_book`). This module adds the *payoff-at-expiry* profile
and break-even solver that traders reason about, plus named constructors for
the common structures:

  * vertical spreads (bull/bear call/put),
  * straddle and strangle,
  * risk reversal (short put, long call — a synthetic-forward skew trade),
  * butterfly and iron condor.

Payoffs are intrinsic values at expiry (per unit, before premium); combined
with the net premium from ``price_book`` they give the full P&L diagram.
"""

import math
from typing import List, Sequence

from .portfolio import Contract, price_book, Book
from .bsm import OptionType, _coerce_type


def _leg(S, K, t, r, sigma, ot, qty, b, mult, label):
    return Contract(S=S, K=K, t=t, r=r, sigma=sigma, option_type=ot,
                    b=b, qty=qty, multiplier=mult, label=label)


def payoff_at_expiry(book: Book, spot_at_expiry: float) -> float:
    """Intrinsic payoff of the book's legs at a terminal spot (per multiplier)."""
    total = 0.0
    for pos in book.positions:
        c = pos.contract
        intrinsic = (max(spot_at_expiry - c.K, 0.0)
                     if c.option_type is OptionType.CALL
                     else max(c.K - spot_at_expiry, 0.0))
        total += c.qty * c.multiplier * intrinsic
    return total


def payoff_profile(book: Book, spots: Sequence[float]) -> List[float]:
    """Payoff at each terminal spot in ``spots``."""
    return [payoff_at_expiry(book, s) for s in spots]


def break_evens(book: Book, lo: float, hi: float, n: int = 2000) -> List[float]:
    """Find terminal spots where total P&L (payoff - net premium) crosses zero.

    Scans ``[lo, hi]`` on a grid and refines each sign change by bisection. Net
    premium is the book's market value now (positive = we paid it).
    """
    premium = book.net.market_value

    def pnl(s):
        return payoff_at_expiry(book, s) - premium

    xs = [lo + (hi - lo) * i / n for i in range(n + 1)]
    roots = []
    prev_s, prev_v = xs[0], pnl(xs[0])
    for s in xs[1:]:
        v = pnl(s)
        if prev_v == 0.0:
            roots.append(prev_s)
        elif prev_v * v < 0:
            a, fa, c_, fc = prev_s, prev_v, s, v
            for _ in range(60):
                m = 0.5 * (a + c_)
                fm = pnl(m)
                if fa * fm <= 0:
                    c_, fc = m, fm
                else:
                    a, fa = m, fm
            roots.append(0.5 * (a + c_))
        prev_s, prev_v = s, v
    return roots


# --- Named strategy constructors ---------------------------------------------
def vertical_spread(S, K_long, K_short, t, r, sigma, kind="call",
                    b=None, mult=1.0):
    """Bull/bear vertical: long one option at K_long, short one at K_short.

    A call spread with K_long < K_short is a bull spread; a put spread with
    K_long > K_short is a bear spread.
    """
    ot = _coerce_type(kind)
    legs = [
        _leg(S, K_long, t, r, sigma, ot, +1, b, mult, f"long {kind} {K_long}"),
        _leg(S, K_short, t, r, sigma, ot, -1, b, mult, f"short {kind} {K_short}"),
    ]
    return price_book(legs)


def straddle(S, K, t, r, sigma, b=None, mult=1.0, qty=1):
    """Long straddle: long a call and a put at the same strike."""
    legs = [
        _leg(S, K, t, r, sigma, OptionType.CALL, qty, b, mult, f"call {K}"),
        _leg(S, K, t, r, sigma, OptionType.PUT, qty, b, mult, f"put {K}"),
    ]
    return price_book(legs)


def strangle(S, K_put, K_call, t, r, sigma, b=None, mult=1.0, qty=1):
    """Long strangle: long an OTM put and an OTM call (K_put < K_call)."""
    legs = [
        _leg(S, K_put, t, r, sigma, OptionType.PUT, qty, b, mult, f"put {K_put}"),
        _leg(S, K_call, t, r, sigma, OptionType.CALL, qty, b, mult, f"call {K_call}"),
    ]
    return price_book(legs)


def risk_reversal(S, K_put, K_call, t, r, sigma, b=None, mult=1.0):
    """Risk reversal: short an OTM put, long an OTM call (a skew/forward trade)."""
    legs = [
        _leg(S, K_put, t, r, sigma, OptionType.PUT, -1, b, mult, f"short put {K_put}"),
        _leg(S, K_call, t, r, sigma, OptionType.CALL, +1, b, mult, f"long call {K_call}"),
    ]
    return price_book(legs)


def butterfly(S, K_low, K_mid, K_high, t, r, sigma, kind="call", b=None, mult=1.0):
    """Long butterfly: +1 K_low, -2 K_mid, +1 K_high (equally spaced strikes)."""
    ot = _coerce_type(kind)
    legs = [
        _leg(S, K_low, t, r, sigma, ot, +1, b, mult, f"{kind} {K_low}"),
        _leg(S, K_mid, t, r, sigma, ot, -2, b, mult, f"{kind} {K_mid} x2"),
        _leg(S, K_high, t, r, sigma, ot, +1, b, mult, f"{kind} {K_high}"),
    ]
    return price_book(legs)


def iron_condor(S, K_put_long, K_put_short, K_call_short, K_call_long,
                t, r, sigma, b=None, mult=1.0):
    """Iron condor: sell an OTM put spread and an OTM call spread.

    Strikes ordered K_put_long < K_put_short < K_call_short < K_call_long.
    Collects premium; profits if the underlying stays between the short strikes.
    """
    legs = [
        _leg(S, K_put_long, t, r, sigma, OptionType.PUT, +1, b, mult, f"long put {K_put_long}"),
        _leg(S, K_put_short, t, r, sigma, OptionType.PUT, -1, b, mult, f"short put {K_put_short}"),
        _leg(S, K_call_short, t, r, sigma, OptionType.CALL, -1, b, mult, f"short call {K_call_short}"),
        _leg(S, K_call_long, t, r, sigma, OptionType.CALL, +1, b, mult, f"long call {K_call_long}"),
    ]
    return price_book(legs)


def ratio_spread(S, K_long, K_short, t, r, sigma, kind="call", ratio=2,
                 b=None, mult=1.0):
    """Ratio spread: long 1 option at K_long, short ``ratio`` at K_short.

    A call ratio spread (K_long < K_short, ratio > 1) is long one lower-strike
    call and short several higher-strike calls -- typically a small credit or
    debit with a capped-profit tent that turns into unlimited downside beyond
    the short strikes. Returns the leg :class:`Book`.
    """
    ot = _coerce_type(kind)
    if ratio < 1:
        raise ValueError("ratio must be >= 1")
    legs = [
        _leg(S, K_long, t, r, sigma, ot, +1, b, mult, f"long {kind} {K_long}"),
        _leg(S, K_short, t, r, sigma, ot, -ratio, b, mult,
             f"short {ratio}x {kind} {K_short}"),
    ]
    return price_book(legs)


def backspread(S, K_short, K_long, t, r, sigma, kind="call", ratio=2,
               b=None, mult=1.0):
    """Backspread: short 1 option at K_short, long ``ratio`` at K_long.

    The mirror of a ratio spread -- net long options, so it profits from a large
    move (unlimited upside for a call backspread) and loses a little in the
    middle. Returns the leg :class:`Book`.
    """
    ot = _coerce_type(kind)
    if ratio < 1:
        raise ValueError("ratio must be >= 1")
    legs = [
        _leg(S, K_short, t, r, sigma, ot, -1, b, mult, f"short {kind} {K_short}"),
        _leg(S, K_long, t, r, sigma, ot, +ratio, b, mult,
             f"long {ratio}x {kind} {K_long}"),
    ]
    return price_book(legs)


def calendar_spread(S, K, t_near, t_far, r, sigma, kind="call",
                    b=None, mult=1.0):
    """Calendar (horizontal) spread: short the near expiry, long the far, same
    strike ``K``.

    A long calendar is short one near-dated option and long one far-dated option
    at the same strike, financed by the faster time decay of the near leg. The
    two legs carry different maturities (``t_near < t_far``), which the
    per-contract ``t`` supports, so ``price_book`` gives the net debit and the
    net Greeks directly. (The expiry payoff diagram is not well defined by
    intrinsics alone, since the far leg still has time value at the near expiry;
    use the net Greeks and price for analysis.) Returns the leg :class:`Book`.
    """
    ot = _coerce_type(kind)
    if not (0.0 < t_near < t_far):
        raise ValueError("require 0 < t_near < t_far")
    legs = [
        _leg(S, K, t_near, r, sigma, ot, -1, b, mult, f"short {kind} {K} @{t_near}"),
        _leg(S, K, t_far, r, sigma, ot, +1, b, mult, f"long {kind} {K} @{t_far}"),
    ]
    return price_book(legs)


def diagonal_spread(S, K_near, K_far, t_near, t_far, r, sigma, kind="call",
                    b=None, mult=1.0):
    """Diagonal spread: short the near expiry at ``K_near``, long the far expiry
    at ``K_far`` -- a calendar with different strikes on the two legs.

    Combines the horizontal (time) and vertical (strike) spreads. Requires
    ``t_near < t_far``; strikes may differ freely. Returns the leg :class:`Book`;
    net price and Greeks come from ``price_book`` (see :func:`calendar_spread` on
    the expiry-payoff caveat).
    """
    ot = _coerce_type(kind)
    if not (0.0 < t_near < t_far):
        raise ValueError("require 0 < t_near < t_far")
    legs = [
        _leg(S, K_near, t_near, r, sigma, ot, -1, b, mult,
             f"short {kind} {K_near} @{t_near}"),
        _leg(S, K_far, t_far, r, sigma, ot, +1, b, mult,
             f"long {kind} {K_far} @{t_far}"),
    ]
    return price_book(legs)


def box_spread(S, K_low, K_high, t, r, sigma, b=None, mult=1.0):
    """Box spread: a bull call spread plus a bear put spread on the same strikes.

    Long call ``K_low`` / short call ``K_high`` (bull call) combined with long
    put ``K_high`` / short put ``K_low`` (bear put). The terminal payoff is the
    constant ``K_high - K_low`` regardless of spot, so the box is a synthetic
    zero-coupon bond: its fair value is the discounted strike width
    ``e^{-rt} (K_high - K_low)``, which the net Greeks confirm are ~zero in spot.
    Returns the leg :class:`Book`.
    """
    if K_high <= K_low:
        raise ValueError("K_high must exceed K_low")
    call = OptionType.CALL
    put = OptionType.PUT
    legs = [
        _leg(S, K_low, t, r, sigma, call, +1, b, mult, f"long call {K_low}"),
        _leg(S, K_high, t, r, sigma, call, -1, b, mult, f"short call {K_high}"),
        _leg(S, K_high, t, r, sigma, put, +1, b, mult, f"long put {K_high}"),
        _leg(S, K_low, t, r, sigma, put, -1, b, mult, f"short put {K_low}"),
    ]
    return price_book(legs)


def synthetic_forward(S, K, t, r, sigma, b=None, mult=1.0):
    """Synthetic long forward: long a call and short a put at the same strike.

    By put-call parity the position replicates a forward struck at ``K``: its
    present value is ``C - P = e^{-bt} S - e^{-rt} K`` (carry ``b``), its delta is
    ~1, and its gamma/vega net to ~zero. Returns the leg :class:`Book`.
    """
    legs = [
        _leg(S, K, t, r, sigma, OptionType.CALL, +1, b, mult, f"long call {K}"),
        _leg(S, K, t, r, sigma, OptionType.PUT, -1, b, mult, f"short put {K}"),
    ]
    return price_book(legs)


def collar(S, K_put, K_call, t, r, sigma, b=None, mult=1.0):
    """Protective collar on a long share: long a put at ``K_put`` (floor) and
    short a call at ``K_call`` (cap), with ``K_put < K_call``.

    Priced here as the two option legs (the underlying share is held
    separately); the net option premium is a small debit or credit depending on
    the skew. The collar caps gains above ``K_call`` and floors losses below
    ``K_put``. Returns the leg :class:`Book`.
    """
    if K_put >= K_call:
        raise ValueError("require K_put < K_call")
    legs = [
        _leg(S, K_put, t, r, sigma, OptionType.PUT, +1, b, mult, f"long put {K_put}"),
        _leg(S, K_call, t, r, sigma, OptionType.CALL, -1, b, mult,
             f"short call {K_call}"),
    ]
    return price_book(legs)


def strategy_report(book, lo=None, hi=None, n=4000):
    """Summarize a strategy's expiry P&L: max profit, max loss, break-evens.

    Scans terminal spots on a grid ``[lo, hi]`` (defaults span a wide range
    around the leg strikes) and returns a dict with the net premium, the maximum
    profit and maximum loss seen on the grid (P&L = payoff - premium), whether
    each is bounded (i.e. not still rising/falling at the grid edge), and the
    break-even spots from :func:`break_evens`.

    P&L is per multiplier, matching :func:`payoff_at_expiry`.
    """
    premium = book.net.market_value
    ks = [pos.contract.K for pos in book.positions]
    S0 = book.positions[0].contract.S if book.positions else 100.0
    if lo is None:
        lo = 0.0
    if hi is None:
        hi = max(ks + [S0]) * 3.0

    xs = [lo + (hi - lo) * i / n for i in range(n + 1)]
    pnls = [payoff_at_expiry(book, s) - premium for s in xs]
    max_profit = max(pnls)
    max_loss = min(pnls)
    # Flag an unbounded tail: still rising/falling at a grid edge and at (near)
    # the extreme there.
    profit_unbounded = ((pnls[-1] >= max_profit - 1e-6 and pnls[-1] > pnls[-2])
                        or (pnls[0] >= max_profit - 1e-6 and pnls[0] > pnls[1]))
    loss_unbounded = ((pnls[0] <= max_loss + 1e-6 and pnls[0] < pnls[1])
                      or (pnls[-1] <= max_loss + 1e-6 and pnls[-1] < pnls[-2]))
    return {
        "net_premium": premium,
        "max_profit": max_profit,
        "max_loss": max_loss,
        "profit_unbounded": profit_unbounded,
        "loss_unbounded": loss_unbounded,
        "break_evens": break_evens(book, lo, hi, n),
    }
