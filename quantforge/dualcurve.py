"""OIS/LIBOR dual-curve discounting.

Post-2008, a floating swap leg is discounted on the collateral (OIS) curve but
its forward rates are *projected* off a separate forward (LIBOR/IBOR) curve --
the two no longer coincide, and their gap is the basis. This module values swaps
under that dual-curve convention:

  * forward rates come from the projection curve,
    ``L(T_{i-1}, T_i) = (P_fwd(T_{i-1})/P_fwd(T_i) - 1) / tau_i``;
  * every cashflow is discounted on the OIS curve;
  * the par swap rate and the value are OIS-discounted, projection-forecast.

It also solves the constant basis spread (added to the projected forwards) that
reprices a set of market par swap rates. Single-curve is the special case where
the two curves are identical -- then the dual-curve par rate collapses to the
single-curve one. Pure standard library, built on :class:`quantforge.DiscountCurve`.
"""

import math

from .discount_curve import DiscountCurve
from .optimize import nelder_mead


def forward_rate(proj_curve, T0, T1):
    """Simply-compounded forward rate over ``[T0, T1]`` off the projection curve."""
    if T1 <= T0:
        raise ValueError("need T1 > T0")
    tau = T1 - T0
    return (proj_curve.df(T0) / proj_curve.df(T1) - 1.0) / tau


def annuity(ois_curve, pay_times):
    """OIS-discounted fixed-leg annuity ``sum_i tau_i P_ois(T_i)``."""
    a = 0.0
    prev = 0.0
    for Ti in pay_times:
        a += (Ti - prev) * ois_curve.df(Ti)
        prev = Ti
    return a


def float_leg_value(ois_curve, proj_curve, pay_times, basis=0.0):
    """Value of a unit-notional float leg: projected forwards, OIS-discounted.

    ``basis`` is a constant spread (in rate units) added to each projected
    forward -- the tenor/currency basis.
    """
    v = 0.0
    prev = 0.0
    for Ti in pay_times:
        tau = Ti - prev
        fwd = forward_rate(proj_curve, prev, Ti) + basis
        v += tau * fwd * ois_curve.df(Ti)
        prev = Ti
    return v


def par_swap_rate(ois_curve, proj_curve, pay_times, basis=0.0):
    """Dual-curve par (fair fixed) rate: float-leg value over the OIS annuity."""
    return float_leg_value(ois_curve, proj_curve, pay_times, basis) \
        / annuity(ois_curve, pay_times)


def swap_value(ois_curve, proj_curve, pay_times, fixed_rate, payer=True,
               basis=0.0):
    """Value of a unit-notional swap (payer = pay fixed, receive float)."""
    flt = float_leg_value(ois_curve, proj_curve, pay_times, basis)
    fixed = fixed_rate * annuity(ois_curve, pay_times)
    return (flt - fixed) if payer else (fixed - flt)


class _ShiftedCurve:
    """A curve whose zero rates are shifted by a constant ``dr`` (continuous).

    ``df_shifted(T) = df(T) * exp(-dr * T)``: a positive ``dr`` raises the zero
    rate at every tenor by ``dr``. Used for parallel-shift risk on a curve that
    exposes only ``df``.
    """

    def __init__(self, curve, dr):
        self._curve = curve
        self._dr = dr

    def df(self, T):
        return self._curve.df(T) * math.exp(-self._dr * T)


def swap_dv01(ois_curve, proj_curve, pay_times, fixed_rate, payer=True,
              basis=0.0, bump=1e-4):
    """Risk of a dual-curve swap by finite differences.

    Returns a dict with:

      * ``pv01`` = the fixed-leg annuity ``sum tau_i P_ois(T_i)`` -- the exact
        ``|dV/d(fixed_rate)|``; a payer's ``dV/d(fixed_rate)`` is ``-pv01``;
      * ``dv01`` = the value change for a 1bp parallel *drop* in **both** curves
        (OIS and projection shifted together), the total delta risk;
      * ``ois_dv01`` / ``proj_dv01`` = the same 1bp-drop risk from shifting only
        the discount (OIS) or only the projection curve.

    ``bump`` is the parallel shift (default 1bp). All figures are per unit
    notional.
    """
    base = swap_value(ois_curve, proj_curve, pay_times, fixed_rate, payer, basis)
    pv01 = annuity(ois_curve, pay_times)

    def val(o, p):
        return swap_value(o, p, pay_times, fixed_rate, payer, basis)

    # 1bp parallel *drop*: zero rates fall by ``bump`` -> dr = -bump.
    ois_dn = _ShiftedCurve(ois_curve, -bump)
    proj_dn = _ShiftedCurve(proj_curve, -bump)
    dv01 = val(ois_dn, proj_dn) - base
    ois_dv01 = val(ois_dn, proj_curve) - base
    proj_dv01 = val(ois_curve, proj_dn) - base
    return {"value": base, "pv01": pv01, "dv01": dv01,
            "ois_dv01": ois_dv01, "proj_dv01": proj_dv01}


def calibrate_basis(ois_curve, proj_curve, swap_maturities, par_rates,
                    freq=1.0):
    """Solve the constant basis spread that reprices the given par swap rates.

    Minimises the squared par-rate error over a single additive basis on the
    projected forwards, with the built-in Nelder-Mead. Returns ``(basis, rmse)``.
    """
    tau = 1.0 / freq
    schedules = []
    for mat in swap_maturities:
        n = int(round(mat * freq))
        schedules.append([(k + 1) * tau for k in range(n)])

    def objective(p):
        basis = p[0]
        err = 0.0
        for pay, target in zip(schedules, par_rates):
            model = par_swap_rate(ois_curve, proj_curve, pay, basis)
            err += (model - target) ** 2
        return err

    best, f = nelder_mead(objective, [0.0], step=0.005, max_iter=2000, tol=1e-18)
    return best[0], math.sqrt(f / len(par_rates))
