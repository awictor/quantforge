"""Equity valuation: cost of capital, dividend discount, and DCF.

Cost of equity via CAPM, the weighted-average cost of capital, the Gordon
constant-growth dividend discount model, a terminal (perpetuity) value, and a
two-stage discounted-cashflow valuation. Pure standard library.
"""


def capm_cost_of_equity(risk_free, beta, market_premium):
    """CAPM cost of equity ``risk_free + beta * market_premium``.

    The required return on equity given its market beta and the equity risk
    premium. Rises with beta.
    """
    return risk_free + beta * market_premium


def wacc(equity_value, debt_value, cost_of_equity, cost_of_debt, tax_rate):
    """Weighted-average cost of capital.

    ``E/V * ke + D/V * kd * (1 - tax)`` with ``V = E + D``. The after-tax blended
    discount rate; lies between the after-tax debt cost and the equity cost.
    """
    v = equity_value + debt_value
    if v <= 0:
        raise ValueError("total capital must be positive")
    if not (0.0 <= tax_rate < 1.0):
        raise ValueError("tax_rate must be in [0, 1)")
    we = equity_value / v
    wd = debt_value / v
    return we * cost_of_equity + wd * cost_of_debt * (1.0 - tax_rate)


def gordon_growth_value(dividend_next, discount_rate, growth):
    """Gordon constant-growth value ``D_1 / (r - g)``.

    The present value of a perpetually-growing dividend. Requires ``r > g``
    (otherwise the sum diverges).
    """
    if discount_rate <= growth:
        raise ValueError("discount_rate must exceed growth")
    return dividend_next / (discount_rate - growth)


def terminal_value(final_cashflow, discount_rate, growth):
    """Gordon terminal (continuing) value at the end of an explicit forecast.

    ``final_cashflow * (1 + growth) / (discount_rate - growth)`` -- the perpetuity
    value of cashflows beyond the forecast horizon. Requires ``r > g``.
    """
    if discount_rate <= growth:
        raise ValueError("discount_rate must exceed growth")
    return final_cashflow * (1.0 + growth) / (discount_rate - growth)


def two_stage_dcf(cashflows, discount_rate, terminal_growth):
    """Two-stage DCF: explicit cashflows plus a discounted terminal value.

    Discounts the explicit ``cashflows`` (periods 1..n) at ``discount_rate`` and
    adds the :func:`terminal_value` of the last cashflow discounted from period
    ``n``. Returns the enterprise/equity value. Requires ``discount_rate >
    terminal_growth``.
    """
    if not cashflows:
        raise ValueError("need at least one cashflow")
    if discount_rate <= terminal_growth:
        raise ValueError("discount_rate must exceed terminal_growth")
    n = len(cashflows)
    pv = sum(cashflows[t] / (1.0 + discount_rate) ** (t + 1) for t in range(n))
    tv = terminal_value(cashflows[-1], discount_rate, terminal_growth)
    pv += tv / (1.0 + discount_rate) ** n
    return pv
