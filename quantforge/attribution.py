"""Greek-based P&L attribution ("P&L explain") for an option position.

Given a position and a market move over one period -- a spot change ``dS``, a
volatility change ``dsigma``, and elapsed time ``dt`` -- this explains the
realized change in the option's value via a second-order Taylor expansion in
the Greeks:

    dP ~= delta*dS + 0.5*gamma*dS^2 + vega*dsigma + theta*dt + rho*dr

The pieces are the *delta*, *gamma*, *vega*, *theta*, and *rho* P&L. The
**unexplained** residual is the true revaluation minus the sum of the pieces --
small for modest moves, larger for big moves (higher-order Greeks). Reporting it
is the point: a large residual flags a move the first/second-order Greeks miss.

Greeks and revaluation both come from the exact BSM engine, so the residual is
a genuine model-consistency check, not an artifact of mismatched pricers.
"""

from dataclasses import dataclass

from .bsm import price, delta, gamma, vega, theta, rho, OptionType, _coerce_type


@dataclass(frozen=True)
class PnLAttribution:
    total: float            # actual revaluation P&L
    delta_pnl: float
    gamma_pnl: float
    vega_pnl: float
    theta_pnl: float
    rho_pnl: float
    explained: float        # sum of the Greek pieces
    unexplained: float      # total - explained


def attribute_pnl(S, K, t, r, sigma, dS, dsigma, dt, dr=0.0,
                  option_type=OptionType.CALL, b=None, qty=1.0) -> PnLAttribution:
    """Attribute an option position's P&L over a move to its Greeks.

    Args:
        dS: change in spot. dsigma: change in vol. dt: elapsed calendar time
            (years). dr: change in rate.
        qty: signed position size (scales every P&L component).

    theta here is the calendar theta (per year) from the engine, so the theta
    P&L is ``theta * dt`` (value lost as time passes). vega is per 1.0 vol, rho
    per 1.0 rate; pass dsigma / dr in those units.
    """
    ot = _coerce_type(option_type)
    if b is None:
        b = r

    # Greeks at the starting point.
    d = delta(S, K, t, r, sigma, ot, b)
    g = gamma(S, K, t, r, sigma, b)
    v = vega(S, K, t, r, sigma, b)
    th = theta(S, K, t, r, sigma, ot, b)
    rh = rho(S, K, t, r, sigma, ot, b)

    delta_pnl = qty * d * dS
    gamma_pnl = qty * 0.5 * g * dS * dS
    vega_pnl = qty * v * dsigma
    theta_pnl = qty * th * dt
    rho_pnl = qty * rh * dr
    explained = delta_pnl + gamma_pnl + vega_pnl + theta_pnl + rho_pnl

    # True revaluation: reprice at the moved market (carry moves with the rate
    # for the plain stock case, matching how rho is defined).
    b_new = b + dr
    p0 = price(S, K, t, r, sigma, ot, b)
    p1 = price(max(S + dS, 1e-12), K, max(t - dt, 1e-12), r + dr,
               max(sigma + dsigma, 1e-12), ot, b_new)
    total = qty * (p1 - p0)

    return PnLAttribution(
        total=total, delta_pnl=delta_pnl, gamma_pnl=gamma_pnl, vega_pnl=vega_pnl,
        theta_pnl=theta_pnl, rho_pnl=rho_pnl, explained=explained,
        unexplained=total - explained,
    )
