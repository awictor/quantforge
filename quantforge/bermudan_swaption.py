"""Bermudan swaption pricing by Longstaff-Schwartz on the G2++ state.

A Bermudan swaption lets the holder enter the underlying swap at any one of a
schedule of exercise dates. Under the two-factor G2++ short-rate model the
state ``(x, y)`` is Markovian and every zero-coupon bond -- hence the swap's
value and the discount factors -- is a closed-form function of that state (see
:mod:`quantforge.g2pp`). So we can:

  1. simulate ``(x, y)`` on the exercise schedule under the risk-neutral measure,
     carrying a discretely-compounded money-market numeraire built from the
     one-period bonds ``P(t_k, t_{k+1})`` along each path;
  2. value the underlying swap at each exercise date in closed form from the
     state; and
  3. run Longstaff-Schwartz backward induction -- regress the discounted
     continuation value on a polynomial basis in ``(x, y)`` and exercise when the
     immediate swap value beats the fitted continuation.

The result is bounded below by the most valuable single European co-terminal
swaption and above by nothing tighter than the full Bermudan, which the tests
check. Pure standard library.
"""

import math
import random

from .g2pp import zero_bond as g2pp_zero_bond, g2pp_V, _B


def _simulate_state(a, b, sigma, eta, rho, times, n_paths, rng):
    """Simulate G2++ (x, y) at the given increasing ``times`` (t_0 = 0 implied).

    Returns ``paths[p]`` = list of ``(x, y)`` at each time in ``times``. Exact
    Gaussian transitions of the OU pair over each step.
    """
    paths = [[] for _ in range(n_paths)]
    prev_t = 0.0
    xs = [0.0] * n_paths
    ys = [0.0] * n_paths
    for tk in times:
        dt = tk - prev_t
        ex = math.exp(-a * dt)
        ey = math.exp(-b * dt)
        vx = sigma * sigma * (1.0 - ex * ex) / (2.0 * a)
        vy = eta * eta * (1.0 - ey * ey) / (2.0 * b)
        # Step covariance of the OU increments.
        cxy = (rho * sigma * eta / (a + b)) * (1.0 - math.exp(-(a + b) * dt))
        sx = math.sqrt(vx)
        sy = math.sqrt(vy)
        corr = cxy / (sx * sy) if sx > 0 and sy > 0 else 0.0
        corr = max(min(corr, 1.0), -1.0)
        for p in range(n_paths):
            z1 = rng.gauss(0.0, 1.0)
            z2 = rng.gauss(0.0, 1.0)
            dxi = sx * z1
            dyi = sy * (corr * z1 + math.sqrt(max(1.0 - corr * corr, 0.0)) * z2)
            xs[p] = xs[p] * ex + dxi
            ys[p] = ys[p] * ey + dyi
            paths[p].append((xs[p], ys[p]))
        prev_t = tk
    return paths


def _swap_value(P0, x, y, a, b, sigma, eta, rho, t, fixed_rate, pay_times, payer):
    """Value at time ``t`` (state x,y) of the swap with fixed leg ``pay_times``.

    ``P0`` is today's discount-factor function ``P(0, .)``. Fixed leg pays
    ``fixed_rate * tau_i`` at each ``pay_times[i]``; the float leg values to
    ``P(t, t) - P(t, T_n)`` at start. A payer swap = float - fixed.
    """
    def Pt(T):
        return g2pp_zero_bond(P0(T), P0(t), x, y, a, b, sigma, eta, rho, t, T)

    start = pay_times[0] - (pay_times[1] - pay_times[0]) if len(pay_times) > 1 \
        else t
    # Assume the swap starts at t (co-terminal exercise): float leg = 1 - P(t,Tn).
    fixed = 0.0
    prev = t
    for Ti in pay_times:
        tau = Ti - prev
        fixed += fixed_rate * tau * Pt(Ti)
        prev = Ti
    float_leg = 1.0 - Pt(pay_times[-1])
    swap = float_leg - fixed
    return swap if payer else -swap


def _poly_regress(X, Y):
    """Least-squares fit of Y on [1, x, y, x^2, xy, y^2]; returns coeff list.

    Solves the 6x6 normal equations by Gaussian elimination (no NumPy).
    """
    def basis(xy):
        x, y = xy
        return [1.0, x, y, x * x, x * y, y * y]
    m = 6
    ATA = [[0.0] * m for _ in range(m)]
    ATy = [0.0] * m
    for xy, yv in zip(X, Y):
        bvec = basis(xy)
        for i in range(m):
            ATy[i] += bvec[i] * yv
            for j in range(m):
                ATA[i][j] += bvec[i] * bvec[j]
    # Gaussian elimination with partial pivoting.
    aug = [ATA[i] + [ATy[i]] for i in range(m)]
    for col in range(m):
        piv = max(range(col, m), key=lambda r: abs(aug[r][col]))
        if abs(aug[piv][col]) < 1e-12:
            aug[col][col] += 1e-9   # ridge for a degenerate column
            piv = col
        aug[col], aug[piv] = aug[piv], aug[col]
        pv = aug[col][col]
        for r in range(m):
            if r == col:
                continue
            f = aug[r][col] / pv
            for c in range(col, m + 1):
                aug[r][c] -= f * aug[col][c]
    coeff = [aug[i][m] / aug[i][i] for i in range(m)]

    def predict(xy):
        return sum(c * bi for c, bi in zip(coeff, basis(xy)))
    return predict


def bermudan_swaption_g2pp(P0, exercise_times, fixed_rate, a, b, sigma, eta, rho,
                           payer=True, n_paths=20000, seed=None):
    """Price a Bermudan swaption under G2++ by Longstaff-Schwartz.

    Args:
        P0: today's discount-factor function ``P(0, T)``.
        exercise_times: increasing list of exercise dates; at date ``t_k`` the
            holder may enter the co-terminal swap paying/receiving ``fixed_rate``
            on the remaining schedule ``exercise_times[k:]`` (annual periods).
        payer: True for a payer swaption (pay fixed), else receiver.
        a, b, sigma, eta, rho: G2++ parameters.

    Returns the Bermudan swaption price (today's value). Uses a
    money-market numeraire built from one-period bonds along each path.
    """
    times = list(exercise_times)
    n_ex = len(times)
    if n_ex < 1:
        raise ValueError("need at least one exercise date")
    rng = random.Random(seed)
    paths = _simulate_state(a, b, sigma, eta, rho, times, n_paths, rng)

    # Money-market numeraire N(t_k) along each path: product of 1/P(t_{j},t_{j+1}).
    # Build discount D(t_k) = 1/N(t_k) from consecutive one-period bonds.
    disc = [[1.0] * n_ex for _ in range(n_paths)]
    for p in range(n_paths):
        d = 1.0
        prev_t = 0.0
        for k, tk in enumerate(times):
            xk, yk = paths[p][k]
            if k == 0:
                # Discount from 0 to t0 uses today's curve bond P(0, t0).
                d = P0(tk)
            else:
                x_prev, y_prev = paths[p][k - 1]
                t_prev = times[k - 1]
                Pstep = g2pp_zero_bond(P0(tk), P0(t_prev), x_prev, y_prev,
                                       a, b, sigma, eta, rho, t_prev, tk)
                d *= Pstep
            disc[p][k] = d
            prev_t = tk

    # Immediate swap value at each exercise date, per path.
    imm = [[0.0] * n_ex for _ in range(n_paths)]
    for k, tk in enumerate(times):
        pay_times = times[k + 1:] + [times[-1] + (times[-1] - times[-2])] \
            if n_ex > 1 else [tk + 1.0]
        # Co-terminal schedule: remaining exercise dates act as fixed-leg pay dates.
        sched = times[k + 1:] if k < n_ex - 1 else [tk + 1.0]
        if not sched:
            sched = [tk + 1.0]
        for p in range(n_paths):
            xk, yk = paths[p][k]
            imm[p][k] = max(_swap_value(P0, xk, yk, a, b, sigma, eta, rho,
                                        tk, fixed_rate, sched, payer), 0.0)

    # LSM backward induction. cashflow[p] = discounted value of the chosen policy.
    cash = [imm[p][n_ex - 1] * disc[p][n_ex - 1] for p in range(n_paths)]
    for k in range(n_ex - 2, -1, -1):
        itm = [p for p in range(n_paths) if imm[p][k] > 0.0]
        if len(itm) >= 8:
            X = [paths[p][k] for p in itm]
            Y = [cash[p] / disc[p][k] for p in itm]   # continuation, undiscounted to t_k
            cont_fn = _poly_regress(X, Y)
            for p in itm:
                cont = cont_fn(paths[p][k])
                if imm[p][k] >= cont:
                    cash[p] = imm[p][k] * disc[p][k]
        # (paths not exercised keep their existing discounted cashflow)
    return sum(cash) / n_paths
