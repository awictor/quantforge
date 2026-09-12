"""Drawdown duration and recovery analytics.

Depth is only half the drawdown story -- how *long* an investor stays underwater
matters as much as how far. From a return series this reports:

  * ``max_drawdown_depth`` -- the deepest peak-to-trough fractional drop;
  * the peak, trough, and recovery indices of that episode;
  * ``longest_underwater`` -- the greatest number of periods spent below a prior
    peak (whether or not the deepest drawdown was also the longest);
  * ``time_to_recovery`` -- periods from the deepest trough back to the prior peak
    (``None`` if it never recovers within the sample).

Complements :func:`~quantforge.max_drawdown` and
:func:`~quantforge.drawdown_curve`. Pure standard library.
"""


def _equity_curve(returns):
    equity = 1.0
    curve = [1.0]
    for r in returns:
        equity *= (1.0 + r)
        curve.append(equity)
    return curve


def drawdown_analytics(returns):
    """Depth-and-duration summary of a return series' drawdowns.

    Returns a dict with ``max_drawdown_depth`` (fractional, non-negative),
    ``peak_index`` / ``trough_index`` / ``recovery_index`` of the deepest episode
    (indices into the equity curve, which has ``len(returns)+1`` points),
    ``time_to_recovery`` (periods from trough to recovery, or ``None`` if never),
    and ``longest_underwater`` (the most periods spent below a prior peak).
    ``recovery_index`` is ``None`` if the deepest drawdown never recovers.
    """
    if not returns:
        return {
            "max_drawdown_depth": 0.0,
            "peak_index": 0,
            "trough_index": 0,
            "recovery_index": 0,
            "time_to_recovery": 0,
            "longest_underwater": 0,
        }

    curve = _equity_curve(returns)
    n = len(curve)

    peak = curve[0]
    peak_idx = 0
    best_depth = 0.0
    best_peak_idx = 0
    best_trough_idx = 0

    # Track the longest underwater stretch (periods strictly below a prior peak).
    longest_uw = 0
    cur_uw = 0

    for i in range(n):
        if curve[i] >= peak:
            peak = curve[i]
            peak_idx = i
            cur_uw = 0
        else:
            cur_uw += 1
            if cur_uw > longest_uw:
                longest_uw = cur_uw
        depth = (peak - curve[i]) / peak
        if depth > best_depth:
            best_depth = depth
            best_peak_idx = peak_idx
            best_trough_idx = i

    # Recovery: first index after the trough that regains the pre-drawdown peak.
    pre_peak_level = curve[best_peak_idx]
    recovery_idx = None
    for i in range(best_trough_idx, n):
        if curve[i] >= pre_peak_level:
            recovery_idx = i
            break

    ttr = None if recovery_idx is None else recovery_idx - best_trough_idx

    return {
        "max_drawdown_depth": best_depth,
        "peak_index": best_peak_idx,
        "trough_index": best_trough_idx,
        "recovery_index": recovery_idx,
        "time_to_recovery": ttr,
        "longest_underwater": longest_uw,
    }
