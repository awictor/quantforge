"""Two-asset PDE solver by the Peaceman-Rachford ADI scheme.

Prices European two-asset options (spread, exchange, basket, ...) by solving the
two-dimensional Black-Scholes PDE in log-prices, where the diffusion
coefficients are constant so the operator splits cleanly. The alternating-
direction-implicit (Peaceman-Rachford) scheme takes each time step in two
half-steps: implicit in ``x1`` then implicit in ``x2``, with the mixed
(correlation) second derivative treated explicitly. Each half-step is a set of
independent tridiagonal solves (Thomas), so the whole step is O(N1 N2).

The grid is uniform in ``x_i = ln S_i`` on a box several standard deviations
around the spots; the far boundaries use the terminal payoff (a reasonable
Dirichlet condition for the smooth two-asset payoffs here). Cross-checked
against the Margrabe exchange closed form. Pure standard library.
"""

import math

from .pde import _thomas


def adi_two_asset(payoff, S1, S2, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0,
                  n1=60, n2=60, n_time=40, width=5.0):
    """Price a European two-asset option by Peaceman-Rachford ADI.

    Args:
        payoff: callable ``payoff(s1, s2)`` giving the terminal value.
        S1, S2: spot prices. sigma1, sigma2, rho: vols and correlation.
        q1, q2: dividend yields. n1, n2, n_time: grid resolutions.
        width: half-width of the log-price box in standard deviations.

    Returns the discounted option value interpolated at ``(S1, S2)``.
    """
    if S1 <= 0 or S2 <= 0:
        raise ValueError("prices must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if t == 0:
        return payoff(S1, S2)

    x1c, x2c = math.log(S1), math.log(S2)
    h1 = width * sigma1 * math.sqrt(t)
    h2 = width * sigma2 * math.sqrt(t)
    x1 = [x1c - h1 + 2.0 * h1 * i / n1 for i in range(n1 + 1)]
    x2 = [x2c - h2 + 2.0 * h2 * j / n2 for j in range(n2 + 1)]
    dx1 = 2.0 * h1 / n1
    dx2 = 2.0 * h2 / n2
    dt = t / n_time

    mu1 = r - q1 - 0.5 * sigma1 * sigma1
    mu2 = r - q2 - 0.5 * sigma2 * sigma2
    a1 = 0.5 * sigma1 * sigma1
    a2 = 0.5 * sigma2 * sigma2
    corr = rho * sigma1 * sigma2

    # Terminal value grid V[i][j].
    V = [[payoff(math.exp(x1[i]), math.exp(x2[j])) for j in range(n2 + 1)]
         for i in range(n1 + 1)]

    def bc(i, j):
        return payoff(math.exp(x1[i]), math.exp(x2[j]))

    # 1D operator coefficients (constant): L1 in x1, L2 in x2, minus 0.5 r each.
    l1_lo = a1 / (dx1 * dx1) - mu1 / (2.0 * dx1)
    l1_di = -2.0 * a1 / (dx1 * dx1) - 0.5 * r
    l1_up = a1 / (dx1 * dx1) + mu1 / (2.0 * dx1)
    l2_lo = a2 / (dx2 * dx2) - mu2 / (2.0 * dx2)
    l2_di = -2.0 * a2 / (dx2 * dx2) - 0.5 * r
    l2_up = a2 / (dx2 * dx2) + mu2 / (2.0 * dx2)

    hdt = 0.5 * dt

    def cross(V, i, j):
        # Explicit mixed second derivative corr * d2V/dx1dx2 (central).
        return corr * (V[i + 1][j + 1] - V[i + 1][j - 1]
                       - V[i - 1][j + 1] + V[i - 1][j - 1]) / (4.0 * dx1 * dx2)

    for _ in range(n_time):
        # ---- Half-step 1: implicit in x1, explicit in x2 + cross. ----
        Vh = [[0.0] * (n2 + 1) for _ in range(n1 + 1)]
        # Boundaries fixed to payoff.
        for i in range(n1 + 1):
            Vh[i][0] = bc(i, 0)
            Vh[i][n2] = bc(i, n2)
        for j in range(n2 + 1):
            Vh[0][j] = bc(0, j)
            Vh[n1][j] = bc(n1, j)
        for j in range(1, n2):
            sub = [0.0] * (n1 + 1)
            dia = [0.0] * (n1 + 1)
            sup = [0.0] * (n1 + 1)
            rhs = [0.0] * (n1 + 1)
            dia[0] = 1.0
            rhs[0] = bc(0, j)
            dia[n1] = 1.0
            rhs[n1] = bc(n1, j)
            for i in range(1, n1):
                # explicit RHS: (I + hdt L2) V + hdt cross.
                expl = (V[i][j]
                        + hdt * (l2_lo * V[i][j - 1] + l2_di * V[i][j]
                                 + l2_up * V[i][j + 1])
                        + hdt * cross(V, i, j))
                sub[i] = -hdt * l1_lo
                dia[i] = 1.0 - hdt * l1_di
                sup[i] = -hdt * l1_up
                rhs[i] = expl
            col = _thomas(sub, dia, sup, rhs)
            for i in range(n1 + 1):
                Vh[i][j] = col[i]

        # ---- Half-step 2: implicit in x2, explicit in x1 + cross. ----
        Vn = [[0.0] * (n2 + 1) for _ in range(n1 + 1)]
        for i in range(n1 + 1):
            Vn[i][0] = bc(i, 0)
            Vn[i][n2] = bc(i, n2)
        for j in range(n2 + 1):
            Vn[0][j] = bc(0, j)
            Vn[n1][j] = bc(n1, j)
        for i in range(1, n1):
            sub = [0.0] * (n2 + 1)
            dia = [0.0] * (n2 + 1)
            sup = [0.0] * (n2 + 1)
            rhs = [0.0] * (n2 + 1)
            dia[0] = 1.0
            rhs[0] = bc(i, 0)
            dia[n2] = 1.0
            rhs[n2] = bc(i, n2)
            for j in range(1, n2):
                expl = (Vh[i][j]
                        + hdt * (l1_lo * Vh[i - 1][j] + l1_di * Vh[i][j]
                                 + l1_up * Vh[i + 1][j])
                        + hdt * cross(Vh, i, j))
                sub[j] = -hdt * l2_lo
                dia[j] = 1.0 - hdt * l2_di
                sup[j] = -hdt * l2_up
                rhs[j] = expl
            row = _thomas(sub, dia, sup, rhs)
            for j in range(n2 + 1):
                Vn[i][j] = row[j]
        V = Vn

    # Bilinear interpolation at (x1c, x2c) -- the centre, exactly a node.
    i = min(max(int((x1c - x1[0]) / dx1), 0), n1 - 1)
    j = min(max(int((x2c - x2[0]) / dx2), 0), n2 - 1)
    wi = (x1c - x1[i]) / dx1
    wj = (x2c - x2[j]) / dx2
    return ((1 - wi) * (1 - wj) * V[i][j] + wi * (1 - wj) * V[i + 1][j]
            + (1 - wi) * wj * V[i][j + 1] + wi * wj * V[i + 1][j + 1])


def adi_spread_option(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0,
                      n1=60, n2=60, n_time=40, width=5.0):
    """European spread-option price ``max(S1 - S2 - K, 0)`` by ADI."""
    disc = math.exp(-r * t)

    def payoff(s1, s2):
        return disc * max(s1 - s2 - K, 0.0)

    # The discount is applied in the payoff; the PDE carries the -0.5 r on each
    # half so e^{-rt} is distributed across the two sweeps. Return raw solve.
    def raw_payoff(s1, s2):
        return max(s1 - s2 - K, 0.0)

    return adi_two_asset(raw_payoff, S1, S2, t, r, sigma1, sigma2, rho,
                         q1, q2, n1, n2, n_time, width)


def adi_two_asset_cs(payoff, S1, S2, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0,
                     n1=60, n2=60, n_time=40, width=5.0, theta=0.5):
    """Price a European two-asset option by the Craig-Sneyd ADI scheme.

    The Peaceman-Rachford scheme (:func:`adi_two_asset`) is only first-order in
    time when a mixed (correlation) derivative is present. Craig-Sneyd fixes that
    with a Douglas predictor followed by a corrector that re-applies the explicit
    cross term at the predicted value, restoring second-order time accuracy. The
    directional operators ``A1``/``A2`` (each carrying half the ``r`` term) are
    solved implicitly (Thomas); the cross term ``A0`` stays explicit.

    Args and return value mirror :func:`adi_two_asset`.
    """
    if S1 <= 0 or S2 <= 0:
        raise ValueError("prices must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if t == 0:
        return payoff(S1, S2)

    x1c, x2c = math.log(S1), math.log(S2)
    h1 = width * sigma1 * math.sqrt(t)
    h2 = width * sigma2 * math.sqrt(t)
    x1 = [x1c - h1 + 2.0 * h1 * i / n1 for i in range(n1 + 1)]
    x2 = [x2c - h2 + 2.0 * h2 * j / n2 for j in range(n2 + 1)]
    dx1 = 2.0 * h1 / n1
    dx2 = 2.0 * h2 / n2
    dt = t / n_time

    mu1 = r - q1 - 0.5 * sigma1 * sigma1
    mu2 = r - q2 - 0.5 * sigma2 * sigma2
    a1 = 0.5 * sigma1 * sigma1
    a2 = 0.5 * sigma2 * sigma2
    corr = rho * sigma1 * sigma2

    # A1 / A2 carry half the discount term each; A0 (cross) carries none.
    l1_lo = a1 / (dx1 * dx1) - mu1 / (2.0 * dx1)
    l1_di = -2.0 * a1 / (dx1 * dx1) - 0.5 * r
    l1_up = a1 / (dx1 * dx1) + mu1 / (2.0 * dx1)
    l2_lo = a2 / (dx2 * dx2) - mu2 / (2.0 * dx2)
    l2_di = -2.0 * a2 / (dx2 * dx2) - 0.5 * r
    l2_up = a2 / (dx2 * dx2) + mu2 / (2.0 * dx2)

    def bc(i, j):
        return payoff(math.exp(x1[i]), math.exp(x2[j]))

    V = [[bc(i, j) for j in range(n2 + 1)] for i in range(n1 + 1)]

    def A1(U, i, j):
        return l1_lo * U[i - 1][j] + l1_di * U[i][j] + l1_up * U[i + 1][j]

    def A2(U, i, j):
        return l2_lo * U[i][j - 1] + l2_di * U[i][j] + l2_up * U[i][j + 1]

    def A0(U, i, j):
        return corr * (U[i + 1][j + 1] - U[i + 1][j - 1]
                       - U[i - 1][j + 1] + U[i - 1][j - 1]) / (4.0 * dx1 * dx2)

    def implicit_x1(Y0, V0):
        """Solve (I - theta dt A1) X = Y0 - theta dt A1 V0, per x2-line."""
        X = [[0.0] * (n2 + 1) for _ in range(n1 + 1)]
        for i in range(n1 + 1):
            X[i][0] = bc(i, 0)
            X[i][n2] = bc(i, n2)
        for j in range(n2 + 1):
            X[0][j] = bc(0, j)
            X[n1][j] = bc(n1, j)
        for j in range(1, n2):
            sub = [0.0] * (n1 + 1)
            dia = [0.0] * (n1 + 1)
            sup = [0.0] * (n1 + 1)
            rhs = [0.0] * (n1 + 1)
            dia[0] = 1.0
            rhs[0] = bc(0, j)
            dia[n1] = 1.0
            rhs[n1] = bc(n1, j)
            for i in range(1, n1):
                sub[i] = -theta * dt * l1_lo
                dia[i] = 1.0 - theta * dt * l1_di
                sup[i] = -theta * dt * l1_up
                rhs[i] = Y0[i][j] - theta * dt * A1(V0, i, j)
            col = _thomas(sub, dia, sup, rhs)
            for i in range(n1 + 1):
                X[i][j] = col[i]
        return X

    def implicit_x2(Y1, V0):
        """Solve (I - theta dt A2) X = Y1 - theta dt A2 V0, per x1-line."""
        X = [[0.0] * (n2 + 1) for _ in range(n1 + 1)]
        for i in range(n1 + 1):
            X[i][0] = bc(i, 0)
            X[i][n2] = bc(i, n2)
        for j in range(n2 + 1):
            X[0][j] = bc(0, j)
            X[n1][j] = bc(n1, j)
        for i in range(1, n1):
            sub = [0.0] * (n2 + 1)
            dia = [0.0] * (n2 + 1)
            sup = [0.0] * (n2 + 1)
            rhs = [0.0] * (n2 + 1)
            dia[0] = 1.0
            rhs[0] = bc(i, 0)
            dia[n2] = 1.0
            rhs[n2] = bc(i, n2)
            for j in range(1, n2):
                sub[j] = -theta * dt * l2_lo
                dia[j] = 1.0 - theta * dt * l2_di
                sup[j] = -theta * dt * l2_up
                rhs[j] = Y1[i][j] - theta * dt * A2(V0, i, j)
            row = _thomas(sub, dia, sup, rhs)
            for j in range(n2 + 1):
                X[i][j] = row[j]
        return X

    for _ in range(n_time):
        # Douglas predictor: explicit full-operator Euler step.
        Y0 = [[V[i][j] for j in range(n2 + 1)] for i in range(n1 + 1)]
        for i in range(1, n1):
            for j in range(1, n2):
                Y0[i][j] = V[i][j] + dt * (A0(V, i, j) + A1(V, i, j)
                                           + A2(V, i, j))
        Y1 = implicit_x1(Y0, V)
        Y2 = implicit_x2(Y1, V)
        # Craig-Sneyd corrector: re-apply the cross term at the predicted value.
        Y0h = [[Y0[i][j] for j in range(n2 + 1)] for i in range(n1 + 1)]
        for i in range(1, n1):
            for j in range(1, n2):
                Y0h[i][j] = Y0[i][j] + 0.5 * dt * (A0(Y2, i, j) - A0(V, i, j))
        Z1 = implicit_x1(Y0h, V)
        Z2 = implicit_x2(Z1, V)
        V = Z2

    i = min(max(int((x1c - x1[0]) / dx1), 0), n1 - 1)
    j = min(max(int((x2c - x2[0]) / dx2), 0), n2 - 1)
    wi = (x1c - x1[i]) / dx1
    wj = (x2c - x2[j]) / dx2
    return ((1 - wi) * (1 - wj) * V[i][j] + wi * (1 - wj) * V[i + 1][j]
            + (1 - wi) * wj * V[i][j + 1] + wi * wj * V[i + 1][j + 1])


def adi_two_asset_american(payoff, S1, S2, t, r, sigma1, sigma2, rho,
                           q1=0.0, q2=0.0, n1=60, n2=60, n_time=40, width=5.0,
                           american=True):
    """Two-asset option with optional early exercise, by Peaceman-Rachford ADI.

    Same discretisation as :func:`adi_two_asset` but, when ``american=True``,
    the value grid is floored at the immediate-exercise payoff after each time
    step (the explicit-payoff projection -- the 2D analogue of the vanilla PSOR
    floor). Suitable for American best-of / worst-of / spread payoffs.

    Returns the value interpolated at ``(S1, S2)``.
    """
    if S1 <= 0 or S2 <= 0:
        raise ValueError("prices must be positive")
    if t < 0:
        raise ValueError("t must be non-negative")
    if t == 0:
        return payoff(S1, S2)

    x1c, x2c = math.log(S1), math.log(S2)
    h1 = width * sigma1 * math.sqrt(t)
    h2 = width * sigma2 * math.sqrt(t)
    x1 = [x1c - h1 + 2.0 * h1 * i / n1 for i in range(n1 + 1)]
    x2 = [x2c - h2 + 2.0 * h2 * j / n2 for j in range(n2 + 1)]
    dx1 = 2.0 * h1 / n1
    dx2 = 2.0 * h2 / n2
    dt = t / n_time

    mu1 = r - q1 - 0.5 * sigma1 * sigma1
    mu2 = r - q2 - 0.5 * sigma2 * sigma2
    a1 = 0.5 * sigma1 * sigma1
    a2 = 0.5 * sigma2 * sigma2
    corr = rho * sigma1 * sigma2

    pay = [[payoff(math.exp(x1[i]), math.exp(x2[j])) for j in range(n2 + 1)]
           for i in range(n1 + 1)]
    V = [row[:] for row in pay]

    l1_lo = a1 / (dx1 * dx1) - mu1 / (2.0 * dx1)
    l1_di = -2.0 * a1 / (dx1 * dx1) - 0.5 * r
    l1_up = a1 / (dx1 * dx1) + mu1 / (2.0 * dx1)
    l2_lo = a2 / (dx2 * dx2) - mu2 / (2.0 * dx2)
    l2_di = -2.0 * a2 / (dx2 * dx2) - 0.5 * r
    l2_up = a2 / (dx2 * dx2) + mu2 / (2.0 * dx2)
    hdt = 0.5 * dt

    def bc(i, j):
        return pay[i][j]

    def cross(U, i, j):
        return corr * (U[i + 1][j + 1] - U[i + 1][j - 1]
                       - U[i - 1][j + 1] + U[i - 1][j - 1]) / (4.0 * dx1 * dx2)

    for _ in range(n_time):
        Vh = [[0.0] * (n2 + 1) for _ in range(n1 + 1)]
        for i in range(n1 + 1):
            Vh[i][0] = bc(i, 0)
            Vh[i][n2] = bc(i, n2)
        for j in range(n2 + 1):
            Vh[0][j] = bc(0, j)
            Vh[n1][j] = bc(n1, j)
        for j in range(1, n2):
            sub = [0.0] * (n1 + 1)
            dia = [0.0] * (n1 + 1)
            sup = [0.0] * (n1 + 1)
            rhs = [0.0] * (n1 + 1)
            dia[0] = 1.0
            rhs[0] = bc(0, j)
            dia[n1] = 1.0
            rhs[n1] = bc(n1, j)
            for i in range(1, n1):
                expl = (V[i][j] + hdt * (l2_lo * V[i][j - 1] + l2_di * V[i][j]
                                         + l2_up * V[i][j + 1])
                        + hdt * cross(V, i, j))
                sub[i] = -hdt * l1_lo
                dia[i] = 1.0 - hdt * l1_di
                sup[i] = -hdt * l1_up
                rhs[i] = expl
            col = _thomas(sub, dia, sup, rhs)
            for i in range(n1 + 1):
                Vh[i][j] = col[i]

        Vn = [[0.0] * (n2 + 1) for _ in range(n1 + 1)]
        for i in range(n1 + 1):
            Vn[i][0] = bc(i, 0)
            Vn[i][n2] = bc(i, n2)
        for j in range(n2 + 1):
            Vn[0][j] = bc(0, j)
            Vn[n1][j] = bc(n1, j)
        for i in range(1, n1):
            sub = [0.0] * (n2 + 1)
            dia = [0.0] * (n2 + 1)
            sup = [0.0] * (n2 + 1)
            rhs = [0.0] * (n2 + 1)
            dia[0] = 1.0
            rhs[0] = bc(i, 0)
            dia[n2] = 1.0
            rhs[n2] = bc(i, n2)
            for j in range(1, n2):
                expl = (Vh[i][j] + hdt * (l1_lo * Vh[i - 1][j] + l1_di * Vh[i][j]
                                          + l1_up * Vh[i + 1][j])
                        + hdt * cross(Vh, i, j))
                sub[j] = -hdt * l2_lo
                dia[j] = 1.0 - hdt * l2_di
                sup[j] = -hdt * l2_up
                rhs[j] = expl
            row = _thomas(sub, dia, sup, rhs)
            for j in range(n2 + 1):
                Vn[i][j] = row[j]
        V = Vn
        if american:
            # Early-exercise projection: floor at the immediate payoff.
            for i in range(n1 + 1):
                for j in range(n2 + 1):
                    if pay[i][j] > V[i][j]:
                        V[i][j] = pay[i][j]

    i = min(max(int((x1c - x1[0]) / dx1), 0), n1 - 1)
    j = min(max(int((x2c - x2[0]) / dx2), 0), n2 - 1)
    wi = (x1c - x1[i]) / dx1
    wj = (x2c - x2[j]) / dx2
    return ((1 - wi) * (1 - wj) * V[i][j] + wi * (1 - wj) * V[i + 1][j]
            + (1 - wi) * wj * V[i][j + 1] + wi * wj * V[i + 1][j + 1])
