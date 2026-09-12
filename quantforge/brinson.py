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
