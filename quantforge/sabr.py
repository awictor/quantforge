"""SABR stochastic-volatility model: Hagan implied-vol expansion + calibration.

The SABR model (Hagan, Kumar, Lesniewski, Woodward, 2002) describes a forward
under stochastic volatility:

    dF = alpha_t * F^beta dW1
    dalpha = nu * alpha_t dW2,   d<W1, W2> = rho dt

with initial vol level ``alpha``, elasticity ``beta`` in [0, 1], vol-of-vol
``nu >= 0`` and correlation ``rho`` in (-1, 1). Hagan's asymptotic formula gives
the Black (lognormal) implied volatility of an option struck at ``K`` on a
forward ``F`` expiring in ``t`` years — the market standard for interpolating
and extrapolating an interest-rate or FX smile.

``sabr_vol`` evaluates the smile; ``calibrate_sabr`` fits (alpha, rho, nu) to a
set of market vols with ``beta`` fixed (the usual convention, since beta and
rho are jointly under-identified from a single smile).
"""

import math
from dataclasses import dataclass
from typing import Sequence, Tuple

from .optimize import nelder_mead


# --------------------------------------------------------------------------
# Forward-mode dual numbers, for exact (machine-precision) parameter
# sensitivities of the Hagan expansion. Each Dual carries a value and a fixed
# vector of partial derivatives w.r.t. the seeded variables
# (F, K, alpha, rho, nu); beta is treated as a constant.
# --------------------------------------------------------------------------
_NVAR = 5  # F, K, alpha, rho, nu


class _Dual:
    __slots__ = ("v", "d")

    def __init__(self, v, d=None):
        self.v = v
        self.d = d if d is not None else [0.0] * _NVAR

    @staticmethod
    def const(v):
        return _Dual(v, [0.0] * _NVAR)

    @staticmethod
    def var(v, i):
        d = [0.0] * _NVAR
        d[i] = 1.0
        return _Dual(v, d)

    def _lift(self, o):
        return o if isinstance(o, _Dual) else _Dual.const(o)

    def __add__(self, o):
        o = self._lift(o)
        return _Dual(self.v + o.v, [a + b for a, b in zip(self.d, o.d)])
    __radd__ = __add__

    def __sub__(self, o):
        o = self._lift(o)
        return _Dual(self.v - o.v, [a - b for a, b in zip(self.d, o.d)])

    def __rsub__(self, o):
        return self._lift(o).__sub__(self)

    def __mul__(self, o):
        o = self._lift(o)
        return _Dual(self.v * o.v,
                     [self.v * b + o.v * a for a, b in zip(self.d, o.d)])
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._lift(o)
        inv = 1.0 / o.v
        return _Dual(self.v * inv,
                     [(a * o.v - self.v * b) * inv * inv
                      for a, b in zip(self.d, o.d)])

    def __rtruediv__(self, o):
        return self._lift(o).__truediv__(self)

    def __pow__(self, p):
        # p is a plain float exponent.
        val = self.v ** p
        coeff = p * self.v ** (p - 1.0)
        return _Dual(val, [coeff * a for a in self.d])

    def log(self):
        inv = 1.0 / self.v
        return _Dual(math.log(self.v), [a * inv for a in self.d])

    def sqrt(self):
        val = math.sqrt(self.v)
        coeff = 0.5 / val
        return _Dual(val, [coeff * a for a in self.d])


def _sabr_vol_dual(F, K, t, alpha, beta, rho, nu):
    """Hagan (2002) vol evaluated on dual numbers; returns the _Dual result.

    Mirrors :func:`sabr_vol` term for term so the value agrees exactly and the
    ``.d`` vector holds the exact partials w.r.t. (F, K, alpha, rho, nu).
    """
    one_beta = 1.0 - beta
    logFK_val = math.log(F.v / K.v)

    if abs(logFK_val) < 1e-12:
        FK_beta = F ** one_beta
        term1 = (one_beta ** 2) / 24.0 * alpha * alpha / (FK_beta ** 2)
        term2 = alpha * (rho * (0.25 * beta * nu)) / FK_beta
        term3 = (nu * nu) * ((2.0 - 3.0 * (rho * rho)) / 24.0)
        return alpha / FK_beta * (1.0 + (term1 + term2 + term3) * t)

    # logFK must stay a dual so the F/K partials propagate through z, x(z) and
    # the log(F/K) prefactor series.
    logFK = (F / K).log()
    FK = F * K
    FK_beta = FK ** (one_beta / 2.0)
    log_FK2 = logFK * logFK

    z = (nu / alpha) * FK_beta * logFK
    inner = (1.0 - 2.0 * rho * z + z * z).sqrt() + z - rho
    x_z = (inner / (1.0 - rho)).log()

    denom = FK_beta * (1.0
                       + (one_beta ** 2) / 24.0 * log_FK2
                       + (one_beta ** 4) / 1920.0 * log_FK2 * log_FK2)

    term1 = (one_beta ** 2) / 24.0 * (alpha * alpha) / (FK ** one_beta)
    term2 = alpha * (rho * (0.25 * beta * nu)) / FK_beta
    term3 = (nu * nu) * ((2.0 - 3.0 * (rho * rho)) / 24.0)
    correction = 1.0 + (term1 + term2 + term3) * t

    return (alpha / denom) * (z / x_z) * correction


def sabr_sensitivities(F, K, t, alpha, beta, rho, nu):
    """Exact partial derivatives of the Hagan SABR vol via forward-mode AD.

    Returns a dict with the vol itself and its machine-precision partials

        ``vol`` and ``d_dF, d_dK, d_dalpha, d_drho, d_dnu``

    computed with dual numbers (no finite-difference truncation error). The
    (alpha, rho, nu) partials are the columns of the calibration Jacobian; the
    ``d_dF`` and ``d_dK`` partials give the smile's backbone and skew slopes.

    At exactly ``F == K`` the ATM branch of :func:`sabr_vol` is used, whose
    F/K partials describe that branch (a finite-difference bump moves off ATM);
    the alpha/rho/nu partials are exact everywhere.
    """
    if F <= 0 or K <= 0:
        raise ValueError("F and K must be positive")
    if alpha <= 0:
        raise ValueError("alpha must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    Fd = _Dual.var(F, 0)
    Kd = _Dual.var(K, 1)
    ad = _Dual.var(alpha, 2)
    rd = _Dual.var(rho, 3)
    nd = _Dual.var(nu, 4)
    out = _sabr_vol_dual(Fd, Kd, t, ad, beta, rd, nd)
    g = out.d
    return {
        "vol": out.v,
        "d_dF": g[0], "d_dK": g[1],
        "d_dalpha": g[2], "d_drho": g[3], "d_dnu": g[4],
    }


def sabr_jacobian(F, t, strikes: Sequence[float], alpha, beta, rho, nu):
    """Calibration Jacobian ``d sabr_vol(K_i) / d (alpha, rho, nu)``.

    Returns a list of ``[d_dalpha, d_drho, d_dnu]`` rows, one per strike, using
    the exact dual-number partials. This is the ``J`` a Gauss-Newton or
    Levenberg-Marquardt step needs, and ``(J^T J)^{-1}`` gives the asymptotic
    parameter covariance for standard errors on a fit.
    """
    rows = []
    for K in strikes:
        s = sabr_sensitivities(F, K, t, alpha, beta, rho, nu)
        rows.append([s["d_dalpha"], s["d_drho"], s["d_dnu"]])
    return rows


@dataclass(frozen=True)
class SABRParams:
    alpha: float
    beta: float
    rho: float
    nu: float


def sabr_vol(F, K, t, alpha, beta, rho, nu) -> float:
    """Hagan (2002) lognormal (Black) implied volatility for the SABR model.

    Uses the standard expansion with the ATM limit handled separately to avoid
    the removable 0/0 singularity at ``F == K``.
    """
    if F <= 0 or K <= 0:
        raise ValueError("F and K must be positive")
    if alpha <= 0:
        raise ValueError("alpha must be positive")
    if t <= 0:
        raise ValueError("t must be positive")

    one_beta = 1.0 - beta
    logFK = math.log(F / K)

    # Common third-order time correction bracket, evaluated at the geometric
    # mean forward-strike; independent of the z/x(z) ratio.
    if abs(logFK) < 1e-12:
        # At-the-money expansion (F == K).
        FK_beta = F ** one_beta
        term1 = (one_beta ** 2) / 24.0 * alpha * alpha / (FK_beta ** 2)
        term2 = 0.25 * rho * beta * nu * alpha / FK_beta
        term3 = (2.0 - 3.0 * rho * rho) / 24.0 * nu * nu
        return alpha / FK_beta * (1.0 + (term1 + term2 + term3) * t)

    FK = F * K
    FK_beta = FK ** (one_beta / 2.0)          # (F K)^{(1-beta)/2}
    log_FK2 = logFK * logFK

    # z and x(z).
    z = (nu / alpha) * FK_beta * logFK
    x_z = math.log((math.sqrt(1.0 - 2.0 * rho * z + z * z) + z - rho) / (1.0 - rho))

    # Prefactor denominator series in log(F/K).
    denom = FK_beta * (1.0
                       + (one_beta ** 2) / 24.0 * log_FK2
                       + (one_beta ** 4) / 1920.0 * log_FK2 * log_FK2)

    term1 = (one_beta ** 2) / 24.0 * alpha * alpha / (FK ** one_beta)
    term2 = 0.25 * rho * beta * nu * alpha / FK_beta
    term3 = (2.0 - 3.0 * rho * rho) / 24.0 * nu * nu
    correction = 1.0 + (term1 + term2 + term3) * t

    return (alpha / denom) * (z / x_z) * correction


def calibrate_sabr(
    F, t, strikes: Sequence[float], market_vols: Sequence[float],
    beta: float = 0.5, weights: Sequence[float] = None,
    initial: SABRParams = None, max_iter: int = 4000,
) -> Tuple[SABRParams, float]:
    """Fit (alpha, rho, nu) of a SABR smile to market Black vols; ``beta`` fixed.

    Returns ``(params, rmse)`` where rmse is the root-mean-square vol error.
    Uses a smooth constrained reparametrization so alpha > 0, nu >= 0 and
    rho in (-1, 1), optimized with the built-in Nelder-Mead.
    """
    strikes = [float(k) for k in strikes]
    market_vols = [float(v) for v in market_vols]
    n = len(strikes)
    if n != len(market_vols) or n < 3:
        raise ValueError("need at least 3 matching (strike, vol) points")
    if weights is None:
        weights = [1.0] * n
    wsum = sum(weights)

    if initial is None:
        # Seed alpha from the ATM-ish vol: alpha ~ vol_atm * F^{1-beta}.
        atm_idx = min(range(n), key=lambda i: abs(strikes[i] - F))
        alpha0 = market_vols[atm_idx] * (F ** (1.0 - beta))
        initial = SABRParams(alpha=max(alpha0, 1e-4), beta=beta, rho=-0.2, nu=0.4)

    def softplus(x):
        return math.log1p(math.exp(-abs(x))) + max(x, 0.0)

    def unpack(p):
        pa, pr, pn = p
        return (softplus(pa) + 1e-8, math.tanh(pr), softplus(pn))

    def inv_softplus(y):
        y = max(y, 1e-8)
        return math.log(math.expm1(y)) if y < 30 else y

    x0 = [inv_softplus(initial.alpha),
          math.atanh(max(min(initial.rho, 0.999), -0.999)),
          inv_softplus(initial.nu)]

    def objective(p):
        alpha, rho, nu = unpack(p)
        err = 0.0
        for i in range(n):
            model = sabr_vol(F, strikes[i], t, alpha, beta, rho, nu)
            diff = model - market_vols[i]
            err += weights[i] * diff * diff
        return err

    best_p, best_f = nelder_mead(objective, x0, step=0.3, max_iter=max_iter, tol=1e-16)
    alpha, rho, nu = unpack(best_p)
    return SABRParams(alpha=alpha, beta=beta, rho=rho, nu=nu), math.sqrt(best_f / wsum)


def _solve3(A, b):
    """Solve a 3x3 linear system ``A x = b`` by Gaussian elimination.

    Returns ``None`` if the matrix is singular; used for the SABR Gauss-Newton
    normal equations where ``A`` is the (damped) 3x3 ``J^T J``.
    """
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(3):
        piv = max(range(col, 3), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-300:
            return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        for r in range(3):
            if r == col:
                continue
            f = M[r][col] / pv
            for c in range(col, 4):
                M[r][c] -= f * M[col][c]
    return [M[i][3] / M[i][i] for i in range(3)]


def calibrate_sabr_lm(
    F, t, strikes: Sequence[float], market_vols: Sequence[float],
    beta: float = 0.5, weights: Sequence[float] = None,
    initial: SABRParams = None, max_iter: int = 100, tol: float = 1e-14,
):
    """Fit (alpha, rho, nu) of a SABR smile by Levenberg-Marquardt.

    Uses the exact analytic Jacobian (:func:`sabr_jacobian`) instead of the
    derivative-free Nelder-Mead of :func:`calibrate_sabr`, so it converges in a
    handful of iterations and lands on the same optimum. Parameters are box
    constrained to the valid region (``alpha > 0``, ``nu >= 0``,
    ``-1 < rho < 1``) by clamping each proposed step.

    Returns ``(params, rmse, n_iter)``.
    """
    strikes = [float(k) for k in strikes]
    market_vols = [float(v) for v in market_vols]
    n = len(strikes)
    if n != len(market_vols) or n < 3:
        raise ValueError("need at least 3 matching (strike, vol) points")
    if weights is None:
        weights = [1.0] * n
    wsum = sum(weights)
    sqrt_w = [math.sqrt(w) for w in weights]

    if initial is None:
        atm_idx = min(range(n), key=lambda i: abs(strikes[i] - F))
        alpha0 = market_vols[atm_idx] * (F ** (1.0 - beta))
        initial = SABRParams(alpha=max(alpha0, 1e-4), beta=beta, rho=-0.2, nu=0.4)

    def clamp(alpha, rho, nu):
        alpha = max(alpha, 1e-8)
        nu = max(nu, 0.0)
        rho = max(min(rho, 0.999), -0.999)
        return alpha, rho, nu

    alpha, rho, nu = clamp(initial.alpha, initial.rho, initial.nu)

    def residuals(alpha, rho, nu):
        return [sqrt_w[i] * (sabr_vol(F, strikes[i], t, alpha, beta, rho, nu)
                             - market_vols[i]) for i in range(n)]

    def sse(res):
        return sum(r * r for r in res)

    res = residuals(alpha, rho, nu)
    cost = sse(res)
    lam = 1e-3
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        # Weighted Jacobian rows d residual_i / d (alpha, rho, nu).
        J = sabr_jacobian(F, t, strikes, alpha, beta, rho, nu)
        J = [[sqrt_w[i] * J[i][k] for k in range(3)] for i in range(n)]

        # Normal-equation pieces: JtJ (3x3) and Jtr (3).
        JtJ = [[sum(J[i][a] * J[i][c] for i in range(n)) for c in range(3)]
               for a in range(3)]
        Jtr = [sum(J[i][a] * res[i] for i in range(n)) for a in range(3)]

        # Levenberg-Marquardt damped step, with lambda adaptation. Grow lambda
        # until a step reduces the cost (toward gradient descent), then accept.
        stepped = False
        for _ in range(30):
            A = [[JtJ[a][c] + (lam if a == c else 0.0) * JtJ[a][a]
                  for c in range(3)] for a in range(3)]
            delta = _solve3(A, [-g for g in Jtr])
            if delta is None:
                lam *= 10.0
                continue
            na, nr, nn = clamp(alpha + delta[0], rho + delta[1], nu + delta[2])
            new_res = residuals(na, nr, nn)
            new_cost = sse(new_res)
            if new_cost < cost:
                alpha, rho, nu = na, nr, nn
                res = new_res
                cost_drop = cost - new_cost
                cost = new_cost
                lam = max(lam * 0.5, 1e-12)
                stepped = True
                break
            lam *= 10.0
        # Stop if no downhill step exists, or the improvement is negligible.
        if not stepped or cost_drop <= tol * (1.0 + cost):
            break

    params = SABRParams(alpha=alpha, beta=beta, rho=rho, nu=nu)
    return params, math.sqrt(cost / wsum), n_iter
