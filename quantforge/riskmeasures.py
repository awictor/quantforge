"""Sample risk measures: VaR, expected shortfall, spectral and entropic risk.

Nonparametric risk measures from a sample of profit-and-loss (or return)
observations. Value-at-risk is the loss quantile; expected shortfall (CVaR)
averages the tail beyond it (a coherent measure); spectral risk weights the whole
loss distribution by a decreasing risk-aversion spectrum; and the entropic risk
measure is the exponential certainty equivalent. Losses are taken as negative
P&L, and each measure is reported as a positive loss magnitude. Pure standard
library.
"""

import math


def _sorted_losses(pnl):
    if not pnl:
        raise ValueError("pnl sample must be non-empty")
    # Losses are negative P&L; sort ascending (worst loss first).
    return sorted(-x for x in pnl)


def value_at_risk(pnl, confidence=0.95):
    """Historical value-at-risk at ``confidence`` (a positive loss magnitude).

    The ``confidence`` quantile of the loss distribution (``-pnl``). Uses the
    lower-index empirical quantile so the VaR is a realized sample loss.
    """
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    losses = _sorted_losses(pnl)
    n = len(losses)
    idx = min(int(math.ceil(confidence * n)) - 1, n - 1)
    idx = max(idx, 0)
    return losses[idx]


def expected_shortfall(pnl, confidence=0.95):
    """Expected shortfall (CVaR): mean loss in the worst ``1 - confidence`` tail.

    Averages the losses at or beyond the VaR quantile. A coherent risk measure,
    always at least the :func:`value_at_risk`.
    """
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be in (0, 1)")
    losses = _sorted_losses(pnl)
    n = len(losses)
    k = max(int(math.ceil((1.0 - confidence) * n)), 1)
    tail = losses[n - k:]
    return sum(tail) / len(tail)


def spectral_risk_exponential(pnl, risk_aversion=5.0):
    """Spectral risk measure with an exponential risk-aversion spectrum.

    Weights the sorted losses by the normalized decreasing spectrum
    ``phi(p) ~ e^{-k(1-p)}`` (heavier weight on worse losses). Coherent for any
    decreasing non-negative spectrum; larger ``risk_aversion`` concentrates weight
    on the tail, raising the measure toward the worst loss.
    """
    if risk_aversion <= 0:
        raise ValueError("risk_aversion must be positive")
    losses = _sorted_losses(pnl)   # ascending: index 0 best-case loss
    n = len(losses)
    # Weight worse losses more: position i (0=best) gets phi at p=(i+1)/n.
    weights = [math.exp(risk_aversion * (i + 1) / n) for i in range(n)]
    total = sum(weights)
    return sum(losses[i] * weights[i] for i in range(n)) / total


def is_subadditive(pnl_a, pnl_b, confidence=0.95):
    """Check the subadditivity ``rho(A+B) <= rho(A) + rho(B)`` for expected shortfall.

    Adds the two P&L series scenario-by-scenario and compares the combined
    expected shortfall to the sum of the standalone ones. Expected shortfall is
    coherent, so this holds (up to a tiny numerical tolerance) for any two aligned
    series -- a diagnostic that diversification never increases ES.
    """
    if len(pnl_a) != len(pnl_b):
        raise ValueError("series must have equal length")
    combined = [a + b for a, b in zip(pnl_a, pnl_b)]
    rho_ab = expected_shortfall(combined, confidence)
    rho_a = expected_shortfall(pnl_a, confidence)
    rho_b = expected_shortfall(pnl_b, confidence)
    return rho_ab <= rho_a + rho_b + 1e-12


def component_expected_shortfall(component_pnls, confidence=0.95):
    """Euler (component) expected-shortfall allocation across sub-portfolios.

    ``component_pnls`` is a list of aligned P&L series whose sum is the portfolio.
    Each component's ES contribution is the average of its own losses over the
    scenarios where the *total* portfolio is in its worst ``1 - confidence`` tail;
    the contributions sum to the portfolio :func:`expected_shortfall` (Euler's
    theorem for the positively-homogeneous ES). Returns the list of contributions.
    """
    if not component_pnls:
        raise ValueError("need at least one component")
    m = len(component_pnls[0])
    if any(len(c) != m for c in component_pnls):
        raise ValueError("all component series must have equal length")
    total = [sum(component_pnls[j][i] for j in range(len(component_pnls)))
             for i in range(m)]
    # Identify the tail scenarios of the total loss.
    losses = [-x for x in total]
    order = sorted(range(m), key=lambda i: losses[i], reverse=True)
    k = max(int(math.ceil((1.0 - confidence) * m)), 1)
    tail = order[:k]
    # Each component's mean loss over those tail scenarios.
    return [sum(-comp[i] for i in tail) / k for comp in component_pnls]


def entropic_risk(pnl, risk_aversion=1.0):
    """Entropic (exponential) risk measure ``(1/theta) ln E[e^{-theta X}]``.

    The exponential-utility certainty equivalent of the loss; ``theta =
    risk_aversion``. Convex and increasing in ``theta``; approaches the mean loss
    ``-E[X]`` as ``theta -> 0`` and the worst loss as ``theta -> inf``.
    """
    if risk_aversion <= 0:
        raise ValueError("risk_aversion must be positive")
    if not pnl:
        raise ValueError("pnl sample must be non-empty")
    n = len(pnl)
    # E[e^{-theta X}] with X the P&L; risk measure on losses is (1/theta) ln E.
    m = max(-risk_aversion * x for x in pnl)   # for numerical stability
    s = sum(math.exp(-risk_aversion * x - m) for x in pnl) / n
    return (m + math.log(s)) / risk_aversion
