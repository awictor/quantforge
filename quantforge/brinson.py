"""Brinson-Hood-Beebower performance attribution.

Decomposes a portfolio's active return (portfolio minus benchmark) versus a
benchmark into per-segment allocation, selection, and interaction effects:

    allocation_i  = (w_p,i - w_b,i) * r_b,i
    selection_i   = w_b,i * (r_p,i - r_b,i)
    interaction_i = (w_p,i - w_b,i) * (r_p,i - r_b,i)

The three effects sum, across segments, to the total active return. Allocation
rewards over/under-weighting segments that beat the benchmark; selection rewards
picking segment holdings that beat their benchmark segment. Pure standard library.
"""


def _validate(pw, bw, pr, br):
    n = len(pw)
    if not (len(bw) == len(pr) == len(br) == n):
        raise ValueError("all segment lists must have equal length")
    if n == 0:
        raise ValueError("need at least one segment")


def allocation_effect(portfolio_weights, benchmark_weights, benchmark_returns):
    """Per-segment allocation effect ``(w_p - w_b) * r_b``.

    Positive when the portfolio overweights segments that outperformed the overall
    benchmark (using the segment benchmark return). Returns one value per segment.
    """
    _validate(portfolio_weights, benchmark_weights, benchmark_returns,
              benchmark_returns)
    return [(pw - bw) * br for pw, bw, br in
            zip(portfolio_weights, benchmark_weights, benchmark_returns)]


def selection_effect(benchmark_weights, portfolio_returns, benchmark_returns):
    """Per-segment selection effect ``w_b * (r_p - r_b)``.

    Positive when the portfolio's holdings within a segment beat that segment's
    benchmark, weighted at the benchmark weight.
    """
    _validate(benchmark_weights, benchmark_weights, portfolio_returns,
              benchmark_returns)
    return [bw * (pr - br) for bw, pr, br in
            zip(benchmark_weights, portfolio_returns, benchmark_returns)]


def interaction_effect(portfolio_weights, benchmark_weights, portfolio_returns,
                       benchmark_returns):
    """Per-segment interaction effect ``(w_p - w_b) * (r_p - r_b)``.

    The cross term -- the combined effect of active weighting and active selection
    in the same segment.
    """
    _validate(portfolio_weights, benchmark_weights, portfolio_returns,
              benchmark_returns)
    return [(pw - bw) * (pr - br) for pw, bw, pr, br in
            zip(portfolio_weights, benchmark_weights, portfolio_returns,
                benchmark_returns)]


def carino_factor(portfolio_return, benchmark_return):
    """Cariño (1999) smoothing factor linking arithmetic effects across periods.

    ``k = ln(1 + r_p) - ln(1 + r_b)) / (r_p - r_b)`` when the returns differ, else
    ``1 / (1 + r_p)``. Scaling each period's arithmetic effects by ``k_t`` and
    dividing by the total-period factor makes the smoothed effects compound
    exactly to the geometric active return -- resolving the residual that plain
    arithmetic summation leaves across multiple periods.
    """
    import math
    if portfolio_return <= -1.0 or benchmark_return <= -1.0:
        raise ValueError("returns must exceed -100%")
    if abs(portfolio_return - benchmark_return) < 1e-12:
        return 1.0 / (1.0 + portfolio_return)
    return (math.log1p(portfolio_return) - math.log1p(benchmark_return)) \
        / (portfolio_return - benchmark_return)


def linked_active_return(portfolio_returns_by_period, benchmark_returns_by_period):
    """Geometrically-linked active return over multiple periods.

    Compounds each side's total return across periods and returns the difference
    of the geometric returns:

        (prod(1 + r_p,t) - 1) - (prod(1 + r_b,t) - 1).

    The quantity multi-period Brinson effects must sum to under Cariño linking.
    """
    n = len(portfolio_returns_by_period)
    if len(benchmark_returns_by_period) != n:
        raise ValueError("period series must have equal length")
    prod_p = 1.0
    prod_b = 1.0
    for rp, rb in zip(portfolio_returns_by_period, benchmark_returns_by_period):
        prod_p *= 1.0 + rp
        prod_b *= 1.0 + rb
    return (prod_p - 1.0) - (prod_b - 1.0)


def carino_linked_effects(period_effects, portfolio_returns_by_period,
                          benchmark_returns_by_period):
    """Cariño-smoothed multi-period effect totals that link geometrically.

    ``period_effects`` is a list of per-period arithmetic effect totals (e.g. the
    allocation totals from :func:`brinson_attribution` each period). Each is scaled
    by its Cariño factor and divided by the total-period factor
    ``k = (ln(1 + R_p) - ln(1 + R_b)) / (R_p - R_b)`` on the compounded returns, so
    the smoothed effects across periods sum to the geometrically-linked active
    return. Returns the smoothed per-period effect list.
    """
    n = len(period_effects)
    if not (len(portfolio_returns_by_period) == len(benchmark_returns_by_period) == n):
        raise ValueError("all period series must have equal length")
    Rp = 1.0
    Rb = 1.0
    for rp, rb in zip(portfolio_returns_by_period, benchmark_returns_by_period):
        Rp *= 1.0 + rp
        Rb *= 1.0 + rb
    Rp -= 1.0
    Rb -= 1.0
    k_total = carino_factor(Rp, Rb)
    out = []
    for eff, rp, rb in zip(period_effects, portfolio_returns_by_period,
                           benchmark_returns_by_period):
        out.append(eff * carino_factor(rp, rb) / k_total)
    return out


def brinson_attribution(portfolio_weights, benchmark_weights, portfolio_returns,
                        benchmark_returns):
    """Full Brinson attribution: allocation, selection, interaction, and totals.

    Returns a dict with per-segment ``allocation``, ``selection``, ``interaction``
    lists, their totals, the ``active_return`` (portfolio minus benchmark total
    return), and ``total_effect`` (allocation + selection + interaction totals).
    The total effect equals the active return by construction.
    """
    alloc = allocation_effect(portfolio_weights, benchmark_weights,
                              benchmark_returns)
    sel = selection_effect(benchmark_weights, portfolio_returns,
                           benchmark_returns)
    inter = interaction_effect(portfolio_weights, benchmark_weights,
                               portfolio_returns, benchmark_returns)
    port_return = sum(w * r for w, r in zip(portfolio_weights, portfolio_returns))
    bench_return = sum(w * r for w, r in zip(benchmark_weights, benchmark_returns))
    total = sum(alloc) + sum(sel) + sum(inter)
    return {
        "allocation": alloc,
        "selection": sel,
        "interaction": inter,
        "allocation_total": sum(alloc),
        "selection_total": sum(sel),
        "interaction_total": sum(inter),
        "active_return": port_return - bench_return,
        "total_effect": total,
    }
