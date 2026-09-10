"""Longstaff-Schwartz least-squares Monte Carlo for Bermudan/American options.

The Longstaff-Schwartz (2001) method prices early-exercise options by
simulation. It steps backward through the exercise dates, and at each date
regresses the discounted future continuation value on a polynomial basis of the
current spot (over the in-the-money paths only). Comparing the fitted
continuation value to the immediate exercise payoff gives the optimal stopping
decision path-by-path.

With enough equally-spaced exercise dates a Bermudan converges to the American
price, which the test suite checks against the binomial tree. Everything is
pure standard library: the polynomial least-squares fit is solved via the
normal equations with a small Gaussian elimination.
"""

import math
import random
from typing import List

from .bsm import OptionType, _coerce_type, _validate


def _solve_normal_equations(X: List[List[float]], y: List[float]) -> List[float]:
    """Least-squares solve of X beta ~= y via the normal equations (X'X)b=X'y.

    ``X`` is n-by-p (rows = samples). Uses Gaussian elimination with partial
    pivoting on the small p-by-p system. Falls back to the mean if singular.
    """
    n = len(X)
    p = len(X[0])
    # Build X'X (p x p) and X'y (p).
    ata = [[0.0] * p for _ in range(p)]
    aty = [0.0] * p
    for i in range(n):
        row = X[i]
        yi = y[i]
        for a in range(p):
            aty[a] += row[a] * yi
            ra = row[a]
            for bcol in range(a, p):
                ata[a][bcol] += ra * row[bcol]
    for a in range(p):
        for bcol in range(a):
            ata[a][bcol] = ata[bcol][a]

    # Gaussian elimination with partial pivoting.
    m = [ata[i][:] + [aty[i]] for i in range(p)]
    for col in range(p):
        piv = max(range(col, p), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-12:
            return [sum(y) / n] + [0.0] * (p - 1)  # singular: constant fit
        m[col], m[piv] = m[piv], m[col]
        pivval = m[col][col]
        for r in range(p):
            if r == col:
                continue
            factor = m[r][col] / pivval
            for c in range(col, p + 1):
                m[r][c] -= factor * m[col][c]
    return [m[i][p] / m[i][i] for i in range(p)]


def bermudan_lsm(S, K, t, r, sigma, option_type=OptionType.PUT, b=None,
                 n_steps=50, n_paths=20_000, degree=3, seed=None) -> float:
    """Price a Bermudan option (exercisable at ``n_steps`` equally-spaced dates).

    Args:
        n_steps: number of exercise opportunities over the life; as it grows the
            price approaches the continuously-exercisable American value.
        degree: polynomial degree of the regression basis in spot.
        b: cost of carry (defaults to r). Dividend yield q enters as b = r - q.

    Returns the option price (in-sample LSM estimate, mildly biased low).
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    rng = random.Random(seed)

    dt = t / n_steps
    drift = (b - 0.5 * sigma * sigma) * dt
    vol = sigma * math.sqrt(dt)
    disc = math.exp(-r * dt)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    # Simulate all paths, storing the spot at each step (columns 0..n_steps).
    paths = [[S] * (n_steps + 1) for _ in range(n_paths)]
    for p in range(n_paths):
        s = S
        for step in range(1, n_steps + 1):
            s *= math.exp(drift + vol * rng.gauss(0.0, 1.0))
            paths[p][step] = s

    def payoff(s):
        return max(sign * (s - K), 0.0)

    # Cashflow at expiry for every path.
    cash = [payoff(paths[p][n_steps]) for p in range(n_paths)]

    # Backward induction over the earlier exercise dates.
    for step in range(n_steps - 1, 0, -1):
        itm = [p for p in range(n_paths) if payoff(paths[p][step]) > 0]
        if len(itm) > degree + 1:
            X = []
            Y = []
            for p in itm:
                s = paths[p][step]
                X.append([s ** d for d in range(degree + 1)])
                Y.append(cash[p] * disc)  # discounted continuation value
            beta = _solve_normal_equations(X, Y)
            for p in itm:
                s = paths[p][step]
                cont = sum(beta[d] * s ** d for d in range(degree + 1))
                exercise = payoff(s)
                if exercise >= cont:
                    cash[p] = exercise
                else:
                    cash[p] = cash[p] * disc
            # Paths not in the money simply discount one step.
            for p in range(n_paths):
                if payoff(paths[p][step]) <= 0:
                    cash[p] *= disc
        else:
            # Too few ITM paths to regress: just discount everything a step.
            for p in range(n_paths):
                cash[p] *= disc

    # Discount the step-1 cashflows back to today.
    return sum(cash) * disc / n_paths


def bermudan_lsm_local_vol(S, K, t, r, local_vol_fn, option_type=OptionType.PUT,
                           q=0.0, n_steps=50, n_paths=20_000, degree=3,
                           seed=None) -> float:
    """Bermudan/American option under a local-volatility surface by LSM.

    Same Longstaff-Schwartz backward induction as :func:`bermudan_lsm`, but each
    Euler step uses the spot- and time-dependent ``local_vol_fn(S, tau)`` (with
    ``tau`` the elapsed forward time) instead of a constant vol -- so it prices
    early-exercise options directly on a calibrated Dupire / SVI local-vol
    surface. Carry is ``b = r - q``. A flat ``local_vol_fn`` reproduces the
    constant-vol LSM price.
    """
    ot = _coerce_type(option_type)
    if S <= 0 or K <= 0:
        raise ValueError("S and K must be positive")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    b = r - q
    rng = random.Random(seed)
    dt = t / n_steps
    sqdt = math.sqrt(dt)
    disc = math.exp(-r * dt)
    sign = 1.0 if ot is OptionType.CALL else -1.0

    paths = [[S] * (n_steps + 1) for _ in range(n_paths)]
    for p in range(n_paths):
        s = S
        tau = 0.0
        for step in range(1, n_steps + 1):
            vol = local_vol_fn(s, tau)
            if vol < 0:
                raise ValueError("local vol must be non-negative")
            s *= math.exp((b - 0.5 * vol * vol) * dt + vol * sqdt * rng.gauss(0.0, 1.0))
            paths[p][step] = s
            tau += dt

    def payoff(s):
        return max(sign * (s - K), 0.0)

    cash = [payoff(paths[p][n_steps]) for p in range(n_paths)]
    for step in range(n_steps - 1, 0, -1):
        itm = [p for p in range(n_paths) if payoff(paths[p][step]) > 0]
        if len(itm) > degree + 1:
            X, Y = [], []
            for p in itm:
                s = paths[p][step]
                X.append([s ** d for d in range(degree + 1)])
                Y.append(cash[p] * disc)
            beta = _solve_normal_equations(X, Y)
            for p in itm:
                s = paths[p][step]
                cont = sum(beta[d] * s ** d for d in range(degree + 1))
                if payoff(s) >= cont:
                    cash[p] = payoff(s)
                else:
                    cash[p] = cash[p] * disc
            for p in range(n_paths):
                if payoff(paths[p][step]) <= 0:
                    cash[p] *= disc
        else:
            for p in range(n_paths):
                cash[p] *= disc
    return sum(cash) * disc / n_paths


def bermudan_max_call_lsm(S1, S2, K, t, r, sigma1, sigma2, rho,
                          q1=0.0, q2=0.0, n_steps=50, n_paths=20_000,
                          seed=None) -> float:
    """American call on the maximum of two assets by Longstaff-Schwartz.

    Prices ``max(max(S1_T, S2_T) - K, 0)`` with early exercise at ``n_steps``
    equally-spaced dates -- the classic two-asset LSM benchmark. Two correlated
    GBMs are simulated (``S2``'s shock is ``rho z1 + sqrt(1-rho^2) z2``) and the
    continuation value is regressed on a quadratic basis in both spots plus the
    running max: ``{1, S1, S2, S1^2, S2^2, S1 S2, max(S1,S2)}``, over the
    in-the-money paths at each date.

    Returns the price (in-sample LSM estimate, mildly biased low). It sits at or
    above the European :func:`quantforge.best_of_call_closed`; with dividends the
    gap is the early-exercise premium.
    """
    ot = OptionType.CALL
    if S1 <= 0 or S2 <= 0 or K <= 0:
        raise ValueError("prices and strike must be positive")
    if t <= 0:
        raise ValueError("t must be positive")
    if not -1.0 <= rho <= 1.0:
        raise ValueError("rho must be in [-1, 1]")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")
    rng = random.Random(seed)

    dt = t / n_steps
    dr1 = (r - q1 - 0.5 * sigma1 * sigma1) * dt
    dr2 = (r - q2 - 0.5 * sigma2 * sigma2) * dt
    v1 = sigma1 * math.sqrt(dt)
    v2 = sigma2 * math.sqrt(dt)
    corr2 = math.sqrt(1.0 - rho * rho)
    disc = math.exp(-r * dt)

    # Simulate both spot paths.
    p1 = [[S1] * (n_steps + 1) for _ in range(n_paths)]
    p2 = [[S2] * (n_steps + 1) for _ in range(n_paths)]
    for p in range(n_paths):
        s1, s2 = S1, S2
        for step in range(1, n_steps + 1):
            z1 = rng.gauss(0.0, 1.0)
            z2 = rng.gauss(0.0, 1.0)
            s1 *= math.exp(dr1 + v1 * z1)
            s2 *= math.exp(dr2 + v2 * (rho * z1 + corr2 * z2))
            p1[p][step] = s1
            p2[p][step] = s2

    def payoff(a, bb):
        return max(max(a, bb) - K, 0.0)

    cash = [payoff(p1[p][n_steps], p2[p][n_steps]) for p in range(n_paths)]

    for step in range(n_steps - 1, 0, -1):
        itm = [p for p in range(n_paths)
               if payoff(p1[p][step], p2[p][step]) > 0]
        if len(itm) > 8:
            X, Y = [], []
            for p in itm:
                a, bb = p1[p][step], p2[p][step]
                X.append([1.0, a, bb, a * a, bb * bb, a * bb, max(a, bb)])
                Y.append(cash[p] * disc)
            beta = _solve_normal_equations(X, Y)
            for p in itm:
                a, bb = p1[p][step], p2[p][step]
                basis = [1.0, a, bb, a * a, bb * bb, a * bb, max(a, bb)]
                cont = sum(beta[j] * basis[j] for j in range(len(basis)))
                ex = payoff(a, bb)
                cash[p] = ex if ex >= cont else cash[p] * disc
            for p in range(n_paths):
                if payoff(p1[p][step], p2[p][step]) <= 0:
                    cash[p] *= disc
        else:
            for p in range(n_paths):
                cash[p] *= disc

    return sum(cash) * disc / n_paths


def bermudan_lsm_greeks(S, K, t, r, sigma, option_type=OptionType.PUT, b=None,
                        n_steps=50, n_paths=40_000, degree=3, seed=None,
                        h_rel=0.01):
    """Delta and gamma of a Bermudan/American LSM price by common-random bumps.

    Prices the option at ``S``, ``S(1 +/- h)`` on the *same* random-number
    stream (each call reseeds ``bermudan_lsm`` with the same ``seed``, so the
    Brownian paths coincide up to the spot scaling and the finite differences
    are low-variance). Returns a dict with ``price``, ``delta`` and ``gamma``
    from central differences; ``h_rel`` is the relative spot bump.

    Common random numbers make the bump estimator far less noisy than
    independent re-pricing; the LSM regression is re-fit at each bump, which is
    the standard practical scheme. Delta is reliable; gamma (a second difference
    over a re-fit regression) is only indicative and needs many paths.
    """
    ot = _coerce_type(option_type)
    _validate(S, K, t, sigma)
    if b is None:
        b = r
    h = h_rel * S
    base = bermudan_lsm(S, K, t, r, sigma, ot, b=b, n_steps=n_steps,
                        n_paths=n_paths, degree=degree, seed=seed)
    up = bermudan_lsm(S + h, K, t, r, sigma, ot, b=b, n_steps=n_steps,
                      n_paths=n_paths, degree=degree, seed=seed)
    dn = bermudan_lsm(S - h, K, t, r, sigma, ot, b=b, n_steps=n_steps,
                      n_paths=n_paths, degree=degree, seed=seed)
    delta = (up - dn) / (2.0 * h)
    gamma = (up - 2.0 * base + dn) / (h * h)
    return {"price": base, "delta": delta, "gamma": gamma}


def bermudan_max_call_lsm_greeks(S1, S2, K, t, r, sigma1, sigma2, rho,
                                 q1=0.0, q2=0.0, n_steps=50, n_paths=40_000,
                                 seed=None, h_rel=0.01):
    """Deltas and cross-gamma of an American max-call by common-random bumps.

    Reprices :func:`bermudan_max_call_lsm` at bumped spots on the *same* seed, so
    the two simulations share their Brownian shocks and the finite differences
    are low-variance. Returns a dict with ``price``, the two spot deltas
    (``delta1`` = dV/dS1, ``delta2`` = dV/dS2), the two own-gammas
    (``gamma1``, ``gamma2``), and the cross-gamma (``cross`` = d2V/dS1 dS2).

    The regression is re-fit at each bump (the standard practical scheme). The
    deltas are reliable; the gammas -- second differences over a re-fit
    regression -- are only indicative and need many paths.
    """
    if seed is None:
        seed = 0
    h1 = h_rel * S1
    h2 = h_rel * S2

    def px(s1, s2):
        return bermudan_max_call_lsm(s1, s2, K, t, r, sigma1, sigma2, rho,
                                     q1, q2, n_steps=n_steps, n_paths=n_paths,
                                     seed=seed)

    base = px(S1, S2)
    u1, d1 = px(S1 + h1, S2), px(S1 - h1, S2)
    u2, d2 = px(S1, S2 + h2), px(S1, S2 - h2)
    delta1 = (u1 - d1) / (2.0 * h1)
    delta2 = (u2 - d2) / (2.0 * h2)
    gamma1 = (u1 - 2.0 * base + d1) / (h1 * h1)
    gamma2 = (u2 - 2.0 * base + d2) / (h2 * h2)
    pp = px(S1 + h1, S2 + h2)
    pm = px(S1 + h1, S2 - h2)
    mp = px(S1 - h1, S2 + h2)
    mm = px(S1 - h1, S2 - h2)
    cross = (pp - pm - mp + mm) / (4.0 * h1 * h2)
    return {"price": base, "delta1": delta1, "delta2": delta2,
            "gamma1": gamma1, "gamma2": gamma2, "cross": cross}
