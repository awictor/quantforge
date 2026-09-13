"""Drawdown-based risk measures: DaR, CDaR, average drawdown.

Standard risk metrics look at return quantiles; drawdown risk looks at the depth of
the peak-to-trough declines an investor actually lives through. From the
drawdown series ``DD_t = 1 - equity_t / running_peak_t``:

- ``average_drawdown`` -- the mean drawdown over the path,
- ``drawdown_at_risk`` -- the ``alpha``-quantile of the drawdown distribution (the
  depth exceeded only ``1 - alpha`` of the time),
- ``conditional_drawdown_at_risk`` -- the mean of the worst ``1 - alpha`` drawdowns
  (Chekhlov-Uryasev-Zabarankin CDaR), a coherent drawdown analogue of expected
  shortfall.

``CDaR >= DaR >= average_drawdown >= 0`` and all are bounded above by the maximum
drawdown. Pure standard library.
"""


def _drawdowns(returns):
    """Fractional drawdown at each step of the compounded equity curve."""
    if not returns:
        raise ValueError("returns must be non-empty")
    equity = 1.0
    peak = 1.0
    dd = []
    for r in returns:
        equity *= (1.0 + r)
        if equity > peak:
            peak = equity
        dd.append(1.0 - equity / peak)
    return dd


def average_drawdown(returns):
    """Mean fractional drawdown over the return path (non-negative)."""
    dd = _drawdowns(returns)
    return sum(dd) / len(dd)


def drawdown_at_risk(returns, confidence=0.95):
    """Drawdown-at-risk: the ``confidence``-quantile of the drawdown distribution.

    The drawdown depth that is exceeded only ``1 - confidence`` of the time. A larger
    confidence gives a deeper (more conservative) threshold. Uses the upper-tail
    order statistic of the drawdown series.
    """
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    dd = sorted(_drawdowns(returns))
    n = len(dd)
    idx = int(confidence * n)
    if idx >= n:
        idx = n - 1
    return dd[idx]


def conditional_drawdown_at_risk(returns, confidence=0.95):
    """Conditional drawdown-at-risk (CDaR): mean of the worst ``1 - confidence`` drawdowns.

    The coherent drawdown analogue of expected shortfall (Chekhlov-Uryasev-
    Zabarankin). Averages the deepest tail of the drawdown distribution beyond the
    :func:`drawdown_at_risk` threshold, so ``CDaR >= DaR``. Falls back to the single
    worst drawdown when the tail holds one observation.
    """
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    dd = sorted(_drawdowns(returns))
    n = len(dd)
    idx = int(confidence * n)
    tail = dd[idx:]
    if not tail:
        tail = [dd[-1]]
    return sum(tail) / len(tail)
