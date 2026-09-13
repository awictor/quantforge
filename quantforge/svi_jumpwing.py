"""Gatheral-Jacquier jump-wing (JW) parameterization of the SVI smile.

The raw SVI parameters ``(a, b, rho, m, s)`` are hard to read. The jump-wing
parameterization (Gatheral & Jacquier, 2014) re-expresses a single expiry ``t`` in
trader-friendly quantities:

- ``v``   : ATM total variance divided by ``t`` (so ``sqrt(v)`` is the ATM vol),
- ``psi`` : the ATM skew ``dw/dk`` at ``k = 0`` scaled by ``1 / (2 sqrt(w_atm))``,
- ``p``   : the slope of the left (put) wing,
- ``c``   : the slope of the right (call) wing,
- ``vtilde`` : the minimum total variance of the slice divided by ``t``.

These map one-to-one to the raw parameters at a fixed ``t``. Pure standard library.
"""

import math
from dataclasses import dataclass

from .svi import SVIParams


@dataclass(frozen=True)
class SVIJumpWing:
    v: float        # ATM variance (per unit time)
    psi: float      # ATM skew
    p: float        # left-wing slope
    c: float        # right-wing slope
    vtilde: float   # minimum variance (per unit time)


def raw_to_jumpwing(params: SVIParams, t: float) -> SVIJumpWing:
    """Convert raw SVI parameters to the jump-wing parameterization at expiry ``t``.

    Uses the closed-form Gatheral-Jacquier map. ``t`` is the expiry in years.
    """
    if t <= 0.0:
        raise ValueError("t must be positive")
    a, b, rho, m, s = params.a, params.b, params.rho, params.m, params.s
    root = math.sqrt(m * m + s * s)
    w_atm = a + b * (-rho * m + root)          # total variance at k = 0
    if w_atm <= 0.0:
        raise ValueError("ATM total variance must be positive")
    v = w_atm / t
    sqrt_w = math.sqrt(w_atm)
    psi = (b / (2.0 * sqrt_w)) * (rho - m / root)
    p = b * (1.0 - rho) / sqrt_w
    c = b * (1.0 + rho) / sqrt_w
    vtilde = (a + b * s * math.sqrt(1.0 - rho * rho)) / t
    return SVIJumpWing(v=v, psi=psi, p=p, c=c, vtilde=vtilde)


def jumpwing_to_raw(jw: SVIJumpWing, t: float) -> SVIParams:
    """Convert jump-wing parameters back to raw SVI at expiry ``t``.

    Inverts :func:`raw_to_jumpwing` with the Gatheral-Jacquier formulas. ``t`` is the
    expiry in years.
    """
    if t <= 0.0:
        raise ValueError("t must be positive")
    v, psi, p, c, vtilde = jw.v, jw.psi, jw.p, jw.c, jw.vtilde
    w = v * t
    if w <= 0.0:
        raise ValueError("ATM variance must be positive")
    sqrt_w = math.sqrt(w)
    b = 0.5 * sqrt_w * (c + p)
    if b < 0.0:
        raise ValueError("wing slopes imply negative b")
    rho = 1.0 - p * sqrt_w / b if b > 0.0 else 0.0
    # From psi = (b / (2 sqrt_w)) (rho - m/root) and the minimum-variance relation.
    beta = rho - 2.0 * psi * sqrt_w / b if b > 0.0 else 0.0
    if abs(beta) > 1.0:
        beta = max(-1.0, min(1.0, beta))
    root_rho = math.sqrt(1.0 - rho * rho)
    if beta == 0.0:
        # No horizontal shift (m = 0); the smile minimum sits at the money.
        m = 0.0
        if b > 0.0 and abs(rho) > 1e-12:
            # v*t = a + b*s and vtilde*t = a + b*s*sqrt(1-rho^2) give s directly.
            s = (v - vtilde) * t / (b * (1.0 - root_rho))
        else:
            # rho = 0 and m = 0: v == vtilde, so the curvature is not recoverable
            # from the jump-wing quantities. Fall back to the ATM variance level.
            s = 0.0
    else:
        alpha = math.copysign(math.sqrt(1.0 / (beta * beta) - 1.0), beta)
        m = (v - vtilde) * t / (b * (-rho + math.copysign(
            math.sqrt(1.0 + alpha * alpha), alpha) - alpha * root_rho))
        s = alpha * m
    if s < 0.0:
        s = 0.0
    a = vtilde * t - b * s * root_rho
    return SVIParams(a=a, b=b, rho=rho, m=m, s=s)
